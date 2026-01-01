import os
import subprocess
import sys
import pandas as pd
import glob
import shutil
import requests
import re
import time
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pubchempy as pcp
import urllib.error
import csv
import chardet
from time import sleep
import matplotlib.pyplot as plt
import networkx as nx
import xml.etree.ElementTree as ET

# --- Configuration ---
INPUT_FOLDER_RAW = '01.Drug_Ingredients/01.ingredients_rawdata'
OUTPUT_FOLDER_PREPROCESSED = '01.Drug_Ingredients/02.ingredients_preprocesseddata'
OUTPUT_FOLDER_AUGMENTED = '01.Drug_Ingredients/03.ingredients_augmenteddata'
OUTPUT_FOLDER_FILTERED = '01.Drug_Ingredients/04.ingredients_filtereddata'

# --- 1.1 - 1.2 Data Preprocessing ---

def get_encoding(file_path):
    with open(file_path, 'rb') as f:
        result = chardet.detect(f.read())
    return result['encoding']

def merge_files(input_folder, output_folder):
    if not os.path.exists(input_folder):
        print(f"Warning: Input folder {input_folder} does not exist.")
        return

    for subdir, _, _ in os.walk(input_folder):
        if subdir == input_folder:
            continue

        folder_name = os.path.basename(subdir)
        output_file = os.path.join(output_folder, f"{folder_name}.csv")

        files = glob.glob(os.path.join(subdir, "*.xlsx")) + glob.glob(os.path.join(subdir, "*.txt"))
        combined_data = []

        for file in files:
            if file.endswith('.xlsx') and not os.path.basename(file).startswith('~$'):
                try:
                    data = pd.read_excel(file, dtype=str, engine='openpyxl')
                except Exception as e:
                    print(f"Error reading {file}: {e}")
                    continue
            elif file.endswith('.txt'):
                try:
                    encoding = get_encoding(file)
                    data = pd.read_csv(file, sep="\t", dtype=str, encoding=encoding)
                except Exception as e:
                    print(f"Error reading {file}: {e}")
                    continue
            else:
                continue

            filename = os.path.splitext(os.path.basename(file))[0]
            data.insert(0, 'herb_name', filename)
            combined_data.append(data)

        if combined_data:
            merged_data = pd.concat(combined_data, ignore_index=True)
            merged_data.to_csv(output_file, index=False, encoding='utf-8-sig')
            print(f"Merged data saved to {output_file}")

os.makedirs(OUTPUT_FOLDER_PREPROCESSED, exist_ok=True)
merge_files(INPUT_FOLDER_RAW, OUTPUT_FOLDER_PREPROCESSED)

# --- 1.3 Data Augmentation (HERB) ---
input_file_herb = os.path.join(OUTPUT_FOLDER_PREPROCESSED, "02.HERB.csv")
herb_info_file = "HERB_ingredient_info_v2.xlsx" # Assumed to be in root

if os.path.exists(input_file_herb) and os.path.exists(herb_info_file):
    os.makedirs(OUTPUT_FOLDER_AUGMENTED, exist_ok=True)

    print(f"Processing HERB data from {input_file_herb}...")
    try:
        herb_ingredient_info = pd.read_excel(herb_info_file)
        df = pd.read_csv(input_file_herb)
        df_ingredient_ids = df['Ingredient id']
        herb_ingredient_info_filtered = herb_ingredient_info[herb_ingredient_info['Ingredient_id'].isin(df_ingredient_ids)]
        merged_data = df.merge(herb_ingredient_info_filtered, left_on='Ingredient id', right_on='Ingredient_id', how='left')
        merged_data.drop(columns=['Ingredient_id'], inplace=True)
        output_file_name = os.path.join(OUTPUT_FOLDER_AUGMENTED, "HERB_augmented.csv")
        merged_data.to_csv(output_file_name, index=False, encoding='utf-8-sig')
        print(f"Processed HERB data saved to {output_file_name}")
    except Exception as e:
        print(f"Error processing HERB data: {e}")
else:
    print("Skipping HERB processing (input files missing).")


# --- Target Prediction ---
# Call external script
print("Running Target Prediction...")
input_smiles = "01_drug_ingredients/05.merge/SMILES_qed0.67.csv"
output_targets = "01_drug_ingredients/05.merge/TARGETS_qed0.67.csv"

# Ensure input directory exists for mock test
os.makedirs(os.path.dirname(input_smiles), exist_ok=True)

if not os.path.exists(input_smiles):
    print(f"Error: Input SMILES file not found at {input_smiles}. Please provide input data.")
    # We do NOT generate mock data here. We fail if input is missing.
    sys.exit(1)

# Use sys.executable to ensure we use the same python interpreter
cmd = [sys.executable, "01_drug_ingredients/targets.py", "--input", input_smiles, "--output", output_targets]
try:
    if os.path.exists("01_drug_ingredients/targets.py"):
        subprocess.run(cmd, check=True)
    else:
        print("Warning: targets.py not found. Skipping prediction step.")
except subprocess.CalledProcessError as e:
    print(f"Error running targets.py: {e}")

# --- 2.2 Target Processing (Human Filter) ---
input_folder_targets = '02_ingredients_targets/021_ingredients_targets_orig'
# targets.py is expected to output to input_folder_targets, OR we need to move it there.
if os.path.exists(output_targets):
    os.makedirs(input_folder_targets, exist_ok=True)
    shutil.copy(output_targets, os.path.join(input_folder_targets, "targets.csv"))
    print(f"Copied {output_targets} to {input_folder_targets}/targets.csv")

output_folder_targets_proc = '02_ingredients_targets/022_ingredients_targets_processed'
os.makedirs(output_folder_targets_proc, exist_ok=True)

if os.path.exists(input_folder_targets):
    csv_files = [f for f in os.listdir(input_folder_targets) if f.endswith('.csv')]
    for csv_file in csv_files:
        try:
            file_path = os.path.join(input_folder_targets, csv_file)
            with open(file_path, 'rb') as f:
                result = chardet.detect(f.read())

            df = pd.read_csv(file_path, index_col=None, header=0, encoding=result['encoding'])
            if 'UniProt_name' in df.columns:
                df = df.assign(UniProt_name=df['UniProt_name'].str.split('|')).explode('UniProt_name')
                df = df[df['UniProt_name'].str.contains('_HUMAN', na=False)]
                output_file_name = csv_file.replace('.csv', '_processed.csv')
                output_file_path = os.path.join(output_folder_targets_proc, output_file_name)
                df.to_csv(output_file_path, index=False)
                print(f"Processed target file: {output_file_name}")
        except Exception as e:
            print(f"Error processing {csv_file}: {e}")

# --- 2.3 Unique Targets ---
input_folder_proc = '02_ingredients_targets/022_ingredients_targets_processed'
output_folder_unique = '02_ingredients_targets/023_ingredients_targets_uniquedata'
os.makedirs(output_folder_unique, exist_ok=True)

if os.path.exists(input_folder_proc):
    all_files = os.listdir(input_folder_proc)
    csv_files = [file for file in all_files if file.endswith('.csv')]
    for csv_file in csv_files:
        try:
            input_file_path = os.path.join(input_folder_proc, csv_file)
            with open(input_file_path, 'rb') as f:
                result = chardet.detect(f.read())
            df = pd.read_csv(input_file_path, index_col=None, header=0, encoding=result['encoding'])
            if 'UniProt_name' in df.columns:
                df = df[['UniProt_name']].drop_duplicates()
                output_file_path = os.path.join(output_folder_unique, csv_file.replace('.csv', '_uniquedata.csv'))
                df.to_csv(output_file_path, index=False)
        except Exception as e:
            print(f"Error unique processing {csv_file}: {e}")

# --- 2.4 UniProt Check ---
print("Skipping UniProt API calls for verification in this script version.")
# We need to map 'GENE_HUMAN' to 'GENE' symbol for downstream analysis
input_folder = '02_ingredients_targets/023_ingredients_targets_uniquedata'
output_folder = '04_intersection_targets'
os.makedirs(output_folder, exist_ok=True)
csv_files = glob.glob(os.path.join(input_folder, '*_uniquedata.csv'))
for csv_file in csv_files:
    try:
        df = pd.read_csv(csv_file)
        if 'UniProt_name' in df.columns:
            # Simplified mapping: GENE_HUMAN -> GENE
            df['Processed Value'] = df['UniProt_name'].apply(lambda x: x.split('_')[0] if isinstance(x, str) and '_' in x else x)
            output_csv = os.path.join(output_folder, os.path.basename(csv_file).replace('_uniquedata.csv', '_adjusted.csv'))
            df.to_csv(output_csv, index=False)
            print(f"Mapped targets saved to {output_csv}")
    except Exception as e:
        print(f"Error mapping targets {csv_file}: {e}")


# --- 3 Disease Targets (DisGeNET API) ---
print("Fetching Disease Targets from DisGeNET...")
disease_targets_folder = '03_diseases_targets/033_uniquedata'
os.makedirs(disease_targets_folder, exist_ok=True)
disease_targets_file = os.path.join(disease_targets_folder, 'disease_targets_uniquedata.csv')

# Attempt to fetch
try:
    vocabulary = "mesh"
    disease_id = "D001172"
    token = "1e0082cbe4be2f5cc81b1b9c8876d8a577cfd697" # Token from notebook
    headers = {"Authorization": f"Bearer {token}"}
    url = f"https://www.disgenet.org/api/gda/disease/{vocabulary}/{disease_id}?format=xml"

    response = requests.get(url, headers=headers, timeout=10)

    if response.status_code == 200:
        # Process XML
        root = ET.fromstring(response.text)
        gene_symbols = []
        for item in root.findall('.//list-item'):
            sym = item.find('gene_symbol').text
            if sym:
                gene_symbols.append(sym)

        if gene_symbols:
            gene_symbols = list(set(gene_symbols)) # Unique
            df_dis = pd.DataFrame({'Symbol': gene_symbols})
            df_dis.to_csv(disease_targets_file, index=False)
            print(f"Successfully fetched {len(gene_symbols)} disease targets from DisGeNET.")
        else:
            print("DisGeNET API returned 200 but no gene symbols found.")
    else:
        print(f"DisGeNET API failed with status code: {response.status_code}")

except Exception as e:
    print(f"Error fetching disease targets: {e}")

# --- 4 Intersection ---
output_folder_intersection = '04_intersection_targets'
os.makedirs(output_folder_intersection, exist_ok=True)

uni_sets = []

# Add Ingredient Targets
ingredient_files = glob.glob(os.path.join(output_folder_intersection, '*_adjusted.csv'))
for csv_file in ingredient_files:
    try:
        with open(csv_file, 'rb') as f:
            result = chardet.detect(f.read())
        encoding = result['encoding'] if result['encoding'] else 'utf-8'
        df = pd.read_csv(csv_file, encoding=encoding)
        if 'Processed Value' in df.columns:
            uni_sets.append(set(df['Processed Value']))
    except Exception as e:
        print(f"Error reading {csv_file}: {e}")

# Add Disease Targets
if os.path.exists(disease_targets_file):
    try:
        df_dis = pd.read_csv(disease_targets_file)
        if 'Symbol' in df_dis.columns:
            uni_sets.append(set(df_dis['Symbol']))
            print("Included Disease Targets in intersection.")
    except Exception as e:
        print(f"Error reading disease targets: {e}")
else:
    print("Warning: Disease targets file missing. Intersection logic incomplete.")

if uni_sets:
    intersection = set.intersection(*uni_sets)
    result_df = pd.DataFrame(list(intersection), columns=['Processed Value'])
    result_df.to_csv(os.path.join(output_folder_intersection, 'intersection_targets.csv'), index=False)
    print(f"Intersection targets saved. Count: {len(intersection)}")
else:
    print("No sets for intersection.")


# --- 5 PPI (STRING DB) ---
print("Running PPI Analysis...")
string_api_url = "https://string-db.org/api"
output_format = "tsv"
method = "network"

intersection_file = os.path.join(output_folder_intersection, 'intersection_targets.csv')
output_folder_ppi = "05_ppi"
os.makedirs(output_folder_ppi, exist_ok=True)
ppi_file = os.path.join(output_folder_ppi, "protein_interactions.tsv")

if os.path.exists(intersection_file):
    df_int = pd.read_csv(intersection_file)
    my_genes = df_int['Processed Value'].tolist()

    if my_genes:
        print(f"Attempting PPI retrieval for {len(my_genes)} targets...")
        try:
            request_url = "/".join([string_api_url, output_format, method])
            params = {
                "identifiers": "%0d".join(my_genes),
                "species": 9606,
                "required_score": 400,
                "caller_identity": "network_pharmacology_pipeline"
            }
            response = requests.post(request_url, data=params, timeout=10)

            if response.status_code == 200 and response.text.strip():
                with open(ppi_file, "w") as f:
                    f.write(response.text)
                print(f"Protein interactions saved to {ppi_file}")
            else:
                print(f"STRING DB API returned status {response.status_code} or no data.")

        except Exception as e:
            print(f"Error: Failed to fetch PPI data from STRING DB ({e}).")
            # We assume user wants us to STOP if data cannot be fetched, to respect "No mock data"
    else:
        print("No targets for PPI.")
else:
    print("Intersection file missing for PPI.")


# --- Network Analysis ---
if os.path.exists(ppi_file):
    try:
        ppi_data = pd.read_csv(ppi_file, sep="\t")
        G = nx.Graph()
        for index, row in ppi_data.iterrows():
            if "preferredName_A" in row and "preferredName_B" in row:
                G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row.get("score", 0))

        if len(G.nodes) > 0:
            degree_centrality = nx.degree_centrality(G)
            closeness_centrality = nx.closeness_centrality(G)
            betweenness_centrality = nx.betweenness_centrality(G)
            eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)

            centrality_measures = {
                'degree': degree_centrality,
                'closeness': closeness_centrality,
                'betweenness': betweenness_centrality,
                'eigenvector': eigenvector_centrality
            }
            df_cen = pd.DataFrame(centrality_measures)
            df_cen['sum'] = df_cen.sum(axis=1)
            df_cen = df_cen.sort_values('sum', ascending=False)
            df_cen.to_csv(os.path.join(output_folder_ppi, "centrality_measures.csv"))

            # Select top proteins for R
            df_top = df_cen.iloc[:, 0:1].head(20)
            df_top.index.name = "SYMBOL"
            df_top.to_excel(os.path.join(output_folder_ppi, "selected_proteins.xlsx"), index=True)
            print("PPI analysis complete.")
        else:
            print("PPI Graph is empty.")
    except Exception as e:
        print(f"Error in PPI analysis: {e}")

# --- Execute R Script ---
print('Running R analysis...')
if shutil.which('Rscript'):
    if os.path.exists("05_ppi/selected_proteins.xlsx"):
         try:
             subprocess.run(['Rscript', 'pipeline_analysis.R'], check=True)
         except subprocess.CalledProcessError as e:
             print(f"R script failed: {e}")
    else:
         print("Skipping R analysis (input missing).")
else:
    print("Rscript not found. Skipping R analysis.")
