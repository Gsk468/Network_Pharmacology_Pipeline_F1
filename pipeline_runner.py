import os
import subprocess
import sys
# Mock get_ipython for !bang commands if we execute blindly, but better to replace them


# --- Cell 0 ---
# Install dependencies - Split to ensure critical packages install even if some optional ones fail
# Core data processing
# Skipped: !pip install pandas numpy openpyxl requests tqdm beautifulsoup4
# Network analysis
# Skipped: !pip install networkx matplotlib
# Target Prediction (Selenium)
# Skipped: !pip install selenium webdriver-manager
# Chemical Informatics
# Skipped: !pip install pubchempy chardet
# Visualization (Optional)
try:
    # Skipped: !pip install matplotlib-venn
    pass
except:
    print('Warning: matplotlib-venn installation failed')
# Bioinformatics (Optional/Advanced) - may fail on newer Python versions due to 'line-profiler' build issues
try:
    # Skipped: !pip install bioservices
    pass
except:
    print('Warning: bioservices installation failed. Some UniProt mapping features might be limited.')


# --- Cell 1 ---
# Import necessary libraries
import os
import glob
import numpy as np
import pandas as pd
import pubchempy as pcp
import time
import urllib.error
import subprocess
import requests
import shutil
import csv
import re
import xml.etree.ElementTree as ET
import openpyxl
import chardet
import json
import sys
from time import sleep
try:
    from matplotlib_venn import venn2
except ImportError:
    print('matplotlib_venn not found, venn diagram plotting might fail.')
import matplotlib.pyplot as plt
import networkx as nx
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
try:
    from bioservices import UniProt
except ImportError:
    print('bioservices not found, uniprot mapping might be limited.')


# --- Cell 2 ---
# Configuration
DISEASE_NAME = "Breast Cancer" # Option to select disease name

# Thresholds for Target Prediction
# Note: SwissTarget and PPB3 use Probability (Higher is better, 0.0-1.0)
#       SEA uses P-Value (Lower is better) and MaxTC (Higher is better)
SWISS_THRESHOLD = 0.0
PPB3_THRESHOLD = 0.0
SEA_PVAL_THRESHOLD = 0.05
SEA_MAXTC_THRESHOLD = 0.57 # Filter for SEA Max Tanimoto Coefficient

print(f"Selected Disease: {DISEASE_NAME}")
print(f"Thresholds: Swiss>={SWISS_THRESHOLD}, PPB3>={PPB3_THRESHOLD}, SEA(P)<={SEA_PVAL_THRESHOLD}, SEA(MaxTC)>={SEA_MAXTC_THRESHOLD}")


# --- Cell 3 ---
# Load data from Sample-Plant-smiles.xlsx and preserve Plant info
import os
import pandas as pd
print('Loading data from Sample-Plant-smiles.xlsx...')
try:
    input_df = pd.read_excel('Sample-Plant-smiles.xlsx')
    print('Data loaded successfully. Head:')
    print(input_df.head())

    # Prepare dataframe for subsequent steps
    processed_df = pd.DataFrame()
    # Map columns: Phytochemical -> Ingredient name, SMILES -> SMILES, CID -> Pubchem_CID
    processed_df['Ingredient name'] = input_df['Phytochemical']
    processed_df['SMILES'] = input_df['SMILES']
    processed_df['Pubchem_CID'] = input_df['CID']
    processed_df['Plant'] = input_df['Plant'] # Preserve Plant column
    processed_df['molecular_formula'] = ''

    # Ensure output directory exists
    input_for_step2 = '01_drug_ingredients/05.merge/no_duplicates/Ingredient_smiles.csv'
    os.makedirs(os.path.dirname(input_for_step2), exist_ok=True)

    processed_df.to_csv(input_for_step2, index=False)
    print(f'Processed data prepared for Step 2: {input_for_step2}')

    # Prepare SMILES list for target_prediction.py
    smiles_file = '01_drug_ingredients/smiles_list.txt'
    with open(smiles_file, 'w') as f:
        for smile in processed_df['SMILES']:
            if pd.notna(smile):
                f.write(str(smile).strip() + '\n')
    print(f'SMILES list prepared for target prediction: {smiles_file}')

except FileNotFoundError:
    print('Error: Sample-Plant-smiles.xlsx not found. Please ensure the file exists.')
except Exception as e:
    print(f'An error occurred: {e}')


# --- Cell 9 ---
input_folder = '02_ingredients_targets/021_ingredients_targets_orig'
output_folder = '02_ingredients_targets/022_ingredients_targets_processed'

# Create output folder (if not exists)
os.makedirs(output_folder, exist_ok=True)

# Get all CSV files in folder

# Ensure input folder exists to avoid FileNotFoundError
if not os.path.exists(input_folder):
    print(f"Warning: Input folder {input_folder} does not exist. Previous step might have failed.")
    os.makedirs(input_folder, exist_ok=True)

csv_files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]

for csv_file in csv_files:
    file_path = os.path.join(input_folder, csv_file)

    # Use chardet to automatically detect encoding format
    with open(file_path, 'rb') as f:
        result = chardet.detect(f.read())

    # Read CSV file
    df = pd.read_csv(file_path, index_col=None, header=0, encoding=result['encoding'])

    # Process UniProt_name column
    df = df.assign(UniProt_name=df['UniProt_name'].str.split('|')).explode('UniProt_name')

    # Only keep rows where UniProt_name column contains *_HUMAN
    df = df[df['UniProt_name'].str.contains('_HUMAN')]

    # Save processed CSV file and rename
    output_file_name = csv_file.replace('.csv', '_processed.csv')
    output_file_path = os.path.join(output_folder, output_file_name)

    df.to_csv(output_file_path, index=False)


# --- Cell 11 ---
input_folder = '02_ingredients_targets/022_ingredients_targets_processed'
output_folder = '02_ingredients_targets/023_ingredients_targets_uniquedata'

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# List all files in the input folder
all_files = os.listdir(input_folder)

# Filter only CSV files
csv_files = [file for file in all_files if file.endswith('.csv')]

for csv_file in csv_files:
    input_file_path = os.path.join(input_folder, csv_file)
    output_file_path = os.path.join(output_folder, csv_file)

    # Detect file encoding
    with open(input_file_path, 'rb') as f:
        result = chardet.detect(f.read())

    # Read file with the detected encoding
    df = pd.read_csv(input_file_path, index_col=None, header=0, encoding=result['encoding'])

    # Drop duplicates in 'UniProt_name' column and keep only unique values
    df = df[['UniProt_name']].drop_duplicates()

    # Save the processed dataframe to a new CSV file
    output_file_path = output_file_path.replace('.csv', '_uniquedata.csv')
    df.to_csv(output_file_path, index=False)


# --- Cell 13 ---
## uniprot rest api species human, status reviewed protein sequences
url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29'
all_fastas = requests.get(url).text


# --- Cell 14 ---
input_folder = '02_ingredients_targets/023_ingredients_targets_uniquedata'
output_folder = '04.Final'
os.makedirs(output_folder, exist_ok=True)

# Read csv file
csv_files = glob.glob(os.path.join(input_folder, '*_uniquedata.csv'))
if csv_files:
    current_file = csv_files[0]
    df = pd.read_csv(current_file)
    # Change export filename
    base_name = os.path.basename(current_file)
    output_csv = os.path.join(output_folder, base_name.replace('_uniquedata.csv', '_adjusted.csv'))
else:
    print('No input files found for Step 2.4'); df = pd.DataFrame({'UniProt_name':[]})
    output_csv = os.path.join(output_folder, 'mock_adjusted.csv')
fasta_list = re.split(r'\n(?=>)', all_fastas)
# Get list of values in UniProt_name column of csv table file
symbol_values = df['UniProt_name'].tolist()
# Build regex pattern
pattern = '|'.join([f'\\b{symbol}\\b' for symbol in symbol_values])
# Match fasta list using regex pattern
matched_fastas = [fasta for fasta in fasta_list if re.search(pattern, fasta)]
# Create DataFrame containing fasta data
fasta_df = pd.DataFrame({'fasta': matched_fastas})

# Except for the header line, only keep content after the second "|" and before the first space in each line
# Process each line of data, extract required content
processed_data = []
original_values = []
for index, row in fasta_df.iterrows():
    value = row['fasta']
    value_parts = value.split('|')
    if len(value_parts) > 2:
        processed_value = value_parts[2].split(' ')[0]
        original_value = [match for match in symbol_values if match in value]
    else:
        processed_value = ''
        original_value = []
    processed_data.append(processed_value)
    original_values.append(original_value)

# Create new DataFrame and store processed data into it
processed_df = pd.DataFrame({
    'Original Value': original_values,
    'Processed Value': processed_data
})

# Export as CSV file
processed_df.to_csv(output_csv, index=False)
print(f"Result exported to {output_csv}")


# --- Cell 16 ---
# Function to search OpenTargets by Disease Name
def fetch_opentargets_data(disease_name, output_file):
    print(f'Searching OpenTargets for: {disease_name}')
    url = 'https://api.platform.opentargets.org/api/v4/graphql'

    # 1. Search for disease ID
    query_search = """
    query Search($queryString: String!) {
      search(queryString: $queryString, entityNames: ["disease"], page: {index: 0, size: 1}) {
        hits {
          id
          name
        }
      }
    }
    """

    try:
        response = requests.post(url, json={"query": query_search, "variables": {"queryString": disease_name}})
        response.raise_for_status()
        data = response.json()
        hits = data.get('data', {}).get('search', {}).get('hits', [])

        if not hits:
            print(f"No disease found for '{disease_name}'")
            return

        disease_id = hits[0]['id']
        print(f"Found disease: {hits[0]['name']} ({disease_id})")

        # 2. Get targets
        query_targets = """
        query diseaseTargets($efoId: String!) {
          disease(efoId: $efoId) {
            associatedTargets(page: {index: 0, size: 200}) { # Limit to 200 for demo
              rows {
                target {
                  approvedSymbol
                }
                score
              }
            }
          }
        }
        """

        response = requests.post(url, json={"query": query_targets, "variables": {"efoId": disease_id}})
        response.raise_for_status()
        data = response.json()

        targets = data.get('data', {}).get('disease', {}).get('associatedTargets', {}).get('rows', [])

        # Format for pipeline: 'symbol', 'overallAssociationScore'
        formatted_data = []
        for row in targets:
            formatted_data.append({
                'symbol': row['target']['approvedSymbol'],
                'overallAssociationScore': row['score']
            })

        df = pd.DataFrame(formatted_data)
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        df.to_csv(output_file, index=False)
        print(f"Saved {len(df)} targets to {output_file}")

    except Exception as e:
        print(f"OpenTargets query failed: {e}")

# Run Search
ot_output = '03_diseases_targets/031_orig/Opentargets_fetched.csv'
fetch_opentargets_data(DISEASE_NAME, ot_output)

# GeneCards Placeholder/Instruction
print(f"For GeneCards, please visit https://www.genecards.org/Search/Keyword?queryString={DISEASE_NAME}")
print("Download the results and save as '03_diseases_targets/031_orig/GeneCards_fetched.csv' if automatic search is not possible.")


# --- Cell 19 ---
# Set export folder path and export filename
output_folder = "03_diseases_targets/031_orig"  # Path can be changed as needed
output_filename = "disgenet_results.xml"  # Filename can be changed as needed
os.makedirs(output_folder, exist_ok=True)
output_file = os.path.join(output_folder, output_filename)

# Input disease information
vocabulary = "mesh" # ICD9CM, ICD10, MeSH, OMIM, DO, EFO, NCI, HPO, MONDO, or ORDO identifier
disease_id = "D001943" # Disease id or list of disease ids separated by "," up to 100.
# Updated Disease ID to D001943 (Breast Cancer) as requested.

# Set API key
token = "1e0082cbe4be2f5cc81b1b9c8876d8a577cfd697" # Please replace with your authorized DisGeNET_REST_API key
headers = {"Authorization": f"Bearer {token}"}

# Send request, get GDA data
url = f"https://www.disgenet.org/api/gda/disease/{vocabulary}/{disease_id}?format=xml"  # Three formats available: TSV, JSON, XML
response = requests.get(url, headers=headers)

# Create folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Output result to file
with open(output_file, "w", encoding="utf-8") as f:
    f.write(response.text)

print(f"Result saved to file: {output_file}")


# --- Cell 22 ---
# Check if directory exists, if not, create it
output_dir = '03_diseases_targets/032_trans'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

try:
    # Parse XML file
    tree = ET.parse('03_diseases_targets/031_orig/disgenet_results.xml')
    root = tree.getroot()

    # Open CSV file and write header row
    with open(os.path.join(output_dir, 'disgenet_results.csv'), 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow(['geneid', 'gene_symbol', 'uniprotid', 'gene_dsi', 'gene_dpi', 'gene_pli', 'protein_class', 'protein_class_name', 'diseaseid', 'disease_name', 'disease_class', 'disease_class_name', 'disease_type', 'disease_semantic_type', 'score', 'ei', 'el', 'year_initial', 'year_final', 'source'])

        # Iterate through XML file and write values to CSV file
        for item in root.findall('.//list-item'):
            geneid = item.find('geneid').text
            gene_symbol = item.find('gene_symbol').text
            uniprotid = item.find('uniprotid').text
            gene_dsi = item.find('gene_dsi').text
            gene_dpi = item.find('gene_dpi').text
            gene_pli = item.find('gene_pli').text
            protein_class = item.find('protein_class').text
            protein_class_name = item.find('protein_class_name').text
            diseaseid = item.find('diseaseid').text
            disease_name = item.find('disease_name').text
            disease_class = item.find('disease_class').text
            disease_class_name = item.find('disease_class_name').text
            disease_type = item.find('disease_type').text
            disease_semantic_type = item.find('disease_semantic_type').text
            score = item.find('score').text
            ei = item.find('ei').text
            el = item.find('el').text
            year_initial = item.find('year_initial').text
            year_final = item.find('year_final').text
            source = item.find('source').text

            # Write values to CSV file
            writer.writerow([geneid, gene_symbol, uniprotid, gene_dsi, gene_dpi, gene_pli, protein_class, protein_class_name, diseaseid, disease_name, disease_class, disease_class_name, disease_type, disease_semantic_type, score, ei, el, year_initial, year_final, source])
except Exception as e:
    print(f'Warning: Failed to process DisGeNET XML: {e}')
    print('Proceeding with other available disease target files (GeneCards, OpenTargets)...')


# --- Cell 24 ---
# Input folder path
input_folder = "03_diseases_targets/031_orig"
output_folder = "03_diseases_targets/032_trans"

# Create output folder
os.makedirs(output_folder, exist_ok=True)

# Iterate through all files
for file in os.listdir(input_folder):
    if file.endswith('.tsv') or file.endswith('.txt'):
        # Read tsv file content
        with open(os.path.join(input_folder, file), 'r', encoding='utf-8') as f:
            tsv_reader = csv.reader(f, delimiter='\t')
            rows = [row for row in tsv_reader]

        # Write tsv file content to csv file
        with open(os.path.join(output_folder, file.replace('.tsv', '.csv').replace('.txt', '.csv')), 'w', newline='', encoding='utf-8') as f:
            csv_writer = csv.writer(f)
            csv_writer.writerows(rows)
    elif file.endswith('.csv'):
        # Directly copy csv file to output folder
        shutil.copyfile(os.path.join(input_folder, file), os.path.join(output_folder, file))
    elif file.endswith('.xlsx'):
        # Read xlsx file
        df = pd.read_excel(os.path.join(input_folder, file))

        # Write xlsx file content to csv file
        csv_file = os.path.join(output_folder, file.replace('.xlsx', '.csv'))
        df.to_csv(csv_file, index=False, encoding='utf-8')


# --- Cell 26 ---
input_folder = '03_diseases_targets/032_trans'
output_folder = '03_diseases_targets/033_screened'

os.makedirs(output_folder, exist_ok=True)

## The following are screening criteria
CTD_inferencescore_threshold = 100 # Screening criteria for files with "CTD" in filename
disgenet_score_threshold = 0.1 # Screening criteria for files with "disgenet" in filename
GeneCards_Relevancescore_threshold = 10 # Screening criteria for files with "GeneCards" in filename
NCBI_org_name = "Homo sapiens" # Screening criteria for files with "NCBI" in filename
Opentargets_overallAssociationScore_threshold = 0.1 # Screening criteria for files with "Opentargets" in filename

for file in os.listdir(input_folder):
    if file.endswith('.csv') and 'CTD' in file:
        output_file = os.path.join(output_folder, file)
        with open(os.path.join(input_folder, file), 'r', newline='', encoding='utf-8') as f_input:
            with open(output_file, 'w', newline='', encoding='utf-8') as f_output:
                reader = csv.reader(f_input)
                writer = csv.writer(f_output)
                header = next(reader)
                writer.writerow(header)

                direct_evidence_index = header.index("Direct Evidence")
                inference_score_index = header.index("Inference Score")

                for row in reader:
                    if row[direct_evidence_index] != "":
                        writer.writerow(row)
                    elif float(row[inference_score_index]) >= CTD_inferencescore_threshold:
                        writer.writerow(row)

    elif file.endswith('.csv') and 'disgenet' in file:
        output_file = os.path.join(output_folder, file)
        with open(os.path.join(input_folder, file), 'r', newline='', encoding='utf-8') as f_input:
            with open(output_file, 'w', newline='', encoding='utf-8') as f_output:
                reader = csv.reader(f_input)
                writer = csv.writer(f_output)
                header = next(reader)
                writer.writerow(header)

                score_index = header.index("score")

                for row in reader:
                    if float(row[score_index]) >= disgenet_score_threshold:
                        writer.writerow(row)

    elif file.endswith('.csv') and 'GeneCards' in file:
        output_file = os.path.join(output_folder, file)
        with open(os.path.join(input_folder, file), 'r', newline='', encoding='utf-8') as f_input:
            with open(output_file, 'w', newline='', encoding='utf-8') as f_output:
                reader = csv.reader(f_input)
                writer = csv.writer(f_output)
                header = next(reader)
                writer.writerow(header)

                relevance_score_index = header.index("Relevance score")

                for row in reader:
                    if float(row[relevance_score_index]) >= GeneCards_Relevancescore_threshold:
                        writer.writerow(row)

    elif file.endswith('.csv') and 'NCBI' in file:
        output_file = os.path.join(output_folder, file)
        with open(os.path.join(input_folder, file), 'r', newline='', encoding='utf-8') as f_input:
            with open(output_file, 'w', newline='', encoding='utf-8') as f_output:
                reader = csv.reader(f_input)
                writer = csv.writer(f_output)
                header = next(reader)
                writer.writerow(header)

                org_name_index = header.index("Org_name")

                for row in reader:
                    if row[org_name_index] == NCBI_org_name:
                        writer.writerow(row)

for file in os.listdir(input_folder):
    if file.endswith('.csv') and 'Opentargets' in file:
        output_file = os.path.join(output_folder, file)
        with open(os.path.join(input_folder, file), 'r', newline='', encoding='utf-8') as f_input:
            with open(output_file, 'w', newline='', encoding='utf-8') as f_output:
                reader = csv.reader(f_input)
                writer = csv.writer(f_output)
                header = next(reader)
                writer.writerow(header)

                overall_score_index = header.index("overallAssociationScore")

                for row in reader:
                    if float(row[overall_score_index]) >= Opentargets_overallAssociationScore_threshold:
                        writer.writerow(row)
    elif file.endswith('.csv') and 'OMIM' in file:
        output_file = os.path.join(output_folder, file)
        shutil.copyfile(os.path.join(input_folder, file), output_file)


# --- Cell 28 ---
# Define input and output folder paths
input_folder = '03_diseases_targets/033_screened'
output_folder = '03_diseases_targets/033_uniquedata'

# Create output folder
os.makedirs(output_folder, exist_ok=True)

# Create an empty list to store all values containing "symbol" column
all_symbol_values = []

# Iterate through all CSV files
for file in os.listdir(input_folder):
    if file.endswith('.csv'):
        # Read CSV file content
        with open(os.path.join(input_folder, file), 'r', encoding='utf-8') as f:
            csv_reader = csv.reader(f)
            header = next(csv_reader)
            symbol_columns = [i for i, col in enumerate(header) if 'symbol' in col.lower()]
            # If "symbol" column exists in CSV file
            if symbol_columns:
                # Iterate through all rows
                for row in csv_reader:
                    # Add values of "symbol" column to all_symbol_values list
                    for col_index in symbol_columns:
                        all_symbol_values.append(row[col_index])

# Remove duplicates from all_symbol_values list, keep only unique values
unique_symbol_values = list(set(all_symbol_values))

# Merge all unique values into a long string
merged_values = ', '.join(unique_symbol_values)

# Write merged values to output file
output_file = os.path.join(output_folder, 'disease_targets_uniquedata.csv')
with open(output_file, 'w', newline='', encoding='utf-8') as f:
    csv_writer = csv.writer(f)
    csv_writer.writerow(['Symbol'])
    for symbol in merged_values.split(','):
        csv_writer.writerow([symbol.strip()])


# --- Cell 31 ---
## uniprot rest api species human, status reviewed protein sequences
url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29'
all_fastas = requests.get(url).text


# --- Cell 33 ---
import re
fasta_list = re.split(r'\n(?=>)', all_fastas)
[fasta for fasta in fasta_list if 'DUS3_HUMAN' in fasta]


# --- Cell 35 ---
input_folder = '02_ingredients_targets/023_ingredients_targets_uniquedata'
output_folder = '04_intersection_targets'
os.makedirs(output_folder, exist_ok=True)

# Read csv file
csv_files = glob.glob(os.path.join(input_folder, '*_uniquedata.csv'))
if csv_files:
    current_file = csv_files[0]
    df = pd.read_csv(current_file)
    # Change export filename
    base_name = os.path.basename(current_file)
    output_csv = os.path.join(output_folder, base_name.replace('_uniquedata.csv', '_adjusted.csv'))
else:
    print('No input files found for Step 2.4'); df = pd.DataFrame({'UniProt_name':[]})
    output_csv = os.path.join(output_folder, 'mock_adjusted.csv')
fasta_list = re.split(r'\n(?=>)', all_fastas)
# Get list of values in UniProt_name column of csv table file
symbol_values = df['UniProt_name'].tolist()
# Build regex pattern
pattern = '|'.join([f'\\b{symbol}\\b' for symbol in symbol_values])
# Match fasta list using regex pattern
matched_fastas = [fasta for fasta in fasta_list if re.search(pattern, fasta)]
# Create DataFrame containing fasta data
fasta_df = pd.DataFrame({'fasta': matched_fastas})

# Except for the header line, only keep content after the second "|" and before the first space in each line
# Process each line of data, extract required content
processed_data = []
original_values = []
for index, row in fasta_df.iterrows():
    value = row['fasta']
    value_parts = value.split('|')
    if len(value_parts) > 2:
        processed_value = value_parts[2].split(' ')[0]
        original_value = [match for match in symbol_values if match in value]
    else:
        processed_value = ''
        original_value = []
    processed_data.append(processed_value)
    original_values.append(original_value)

# Create new DataFrame and store processed data into it
processed_df = pd.DataFrame({
    'Original Value': original_values,
    'Processed Value': processed_data
})

# Export as CSV file
processed_df.to_csv(output_csv, index=False)
print(f"Result exported to {output_csv}")


# --- Cell 37 ---
csv_file = '03_diseases_targets/033_uniquedata/disease_targets_uniquedata.csv'
output_csv = '04_intersection_targets/disease_targets_adjusted.csv'

# Match entry name of fasta file by gene name
# Read csv file
df = pd.read_csv(csv_file)
fasta_list = re.split(r'\n(?=>)', all_fastas)
# Get list of values in symbol column of csv table file
symbol_values = df['Symbol'].tolist()
# Build regex pattern
pattern = '|'.join([f'GN={symbol} ' for symbol in symbol_values])

# Match fasta list using regex pattern
matched_fastas = []
matched_symbols = []
for fasta in fasta_list:
    match = re.search(pattern, fasta)
    if match:
        symbol = match.group().split('=')[1]  # Get matched gene symbol
        matched_fastas.append(fasta)
        matched_symbols.append(symbol)

# Create DataFrame containing fasta data
fasta_df = pd.DataFrame({'fasta': matched_fastas, 'Symbol': matched_symbols})

# Except for the header line, only keep content after the second "|" and before the first space in each line
# Process each line of data, extract required content
processed_data = []
for index, row in fasta_df.iterrows():
    value = row['fasta']
    value_parts = value.split('|')
    if len(value_parts) > 2:
        processed_value = value_parts[2].split(' ')[0]
    else:
        processed_value = ''
    processed_data.append(processed_value)

# Create new DataFrame and store processed data and matched gene symbol into it
processed_df = pd.DataFrame({
    'Original Gene Symbol': fasta_df['Symbol'],
    'Processed Value': processed_data
})

# Export as CSV file
processed_df.to_csv(output_csv, index=False)
print(f"Result exported to {output_csv}")


# --- Cell 39 ---
def detect_encoding(file_path):
    with open(file_path, 'rb') as f:
        result = chardet.detect(f.read())
    return result['encoding']

path = '04_intersection_targets'
output_path = '04_intersection_targets'
csv_files = glob.glob(os.path.join(path, '*.csv'))

uni_sets = []

for csv_file in csv_files:
    encoding = detect_encoding(csv_file)
    df = pd.read_csv(csv_file, encoding=encoding)
    uni_set = set(df['Processed Value'])
    uni_sets.append(uni_set)

intersection = set.intersection(*uni_sets)
result_df = pd.DataFrame(list(intersection), columns=['Processed Value'])
result_df.to_csv(os.path.join(output_path, 'intersection_targets.csv'), index=False)

# Plot Venn diagram (only applicable for 2 CSV files)
if len(uni_sets) == 2:
    venn = venn2([uni_sets[0], uni_sets[1]], set_labels=['File 1', 'File 2'])
    plt.savefig(os.path.join(output_path, 'venn_diagram.svg'), format='svg')
    plt.show()
else:
    print("Venn diagram can only be plotted for 2 CSV files.")


# --- Cell 43 ---
output_folder = "05_ppi" # Define output folder path
os.makedirs(output_folder, exist_ok=True)

string_api_url = "https://version-11-5.string-db.org/api"
output_format = "svg"
method = "network"

df = pd.read_csv("04_intersection_targets/intersection_targets.csv") # Read CSV file
my_genes = df["Processed Value"].tolist()

identifiers = "%0d".join(my_genes)

request_url = "/".join([string_api_url, output_format, method])

params = {
    "identifiers": identifiers,
    "species": 9606,
    "required_score": 900,  # "Minimum required interaction score" in STRING webpage, optional 400, 700, 900
    "hide_disconnected_nodes": 1,
    "caller_identity": "www.my_app.org"
}

response = requests.post(request_url, data=params)

file_name = "all_genes_network.svg"
output_path = os.path.join(output_folder, file_name)
print(f"Saving interaction network to {output_path}")

with open(output_path, 'wb') as fh:
    fh.write(response.content)

sleep(1)


# --- Cell 45 ---
output_folder = "05_ppi" # Please modify output folder path according to actual situation
os.makedirs(output_folder, exist_ok=True)

string_api_url = "https://string-db.org/api"
output_format = "tsv"
method = "network"

csv_file_path = "04_intersection_targets/intersection_targets.csv"  # Please modify CSV file path according to actual situation
df = pd.read_csv(csv_file_path)
my_genes = df["Processed Value"].tolist()

request_url = "/".join([string_api_url, output_format, method])

params = {
    "identifiers": "%0d".join(my_genes),
    "species": 9606,
    "required_score": 900, # "Minimum required interaction score" in STRING webpage, optional 400, 700, 900
    "caller_identity": "www.my_example_app.org",
}

response = requests.post(request_url, data=params)


# Save as TSV file
output_file_path = os.path.join(output_folder, "protein_interactions.tsv")
with open(output_file_path, "w") as f:
    f.write(response.text)

print(f"Protein interactions saved to {output_file_path}")


# --- Cell 47 ---
output_folder = "05_ppi"
os.makedirs(output_folder, exist_ok=True)
# Create an empty undirected weighted network
G = nx.Graph()

# 1. Add PPI edges (Target-Target)
ppi_file = "05_ppi/protein_interactions.tsv"
if os.path.exists(ppi_file):
    try:
        ppi_data = pd.read_csv(ppi_file, sep="\t")
        if "preferredName_A" in ppi_data.columns and "preferredName_B" in ppi_data.columns:
            for index, row in ppi_data.iterrows():
                G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row.get("score", 0.0), type='ppi')
        else:
            print("PPI file missing required columns. Skipping PPI edges.")
    except Exception as e:
        print(f"Error reading PPI file: {e}")
else:
    print("PPI file not found. Skipping PPI edges.")

# 2. Add Ingredient-Target edges
# Load targets mapping
try:
    ing_target_file = "02_ingredients_targets/021_ingredients_targets_orig/targets.csv"
    if os.path.exists(ing_target_file):
        it_df = pd.read_csv(ing_target_file)

        # We need to intersect with valid targets from Step 4
        intersection_file = "04_intersection_targets/intersection_targets.csv"
        valid_targets = set()
        if os.path.exists(intersection_file):
             idf = pd.read_csv(intersection_file)
             if 'Processed Value' in idf.columns:
                 valid_targets = set(idf['Processed Value'])

        for idx, row in it_df.iterrows():
            ing = row['Ingredient name']
            targs = str(row['UniProt_name']).split('|')
            for t in targs:
                # Add edge if target is in intersection (or if we want full network)
                # Generally we only care about intersection targets
                if t in valid_targets:
                     G.add_edge(ing, t, weight=1.0, type='ing-target')
                elif not valid_targets: # If intersection empty/missing, add all (fallback)
                     G.add_edge(ing, t, weight=1.0, type='ing-target')

    # 3. Add Plant-Ingredient edges
    ing_smiles_file = "01_drug_ingredients/05.merge/no_duplicates/Ingredient_smiles.csv"
    if os.path.exists(ing_smiles_file):
        is_df = pd.read_csv(ing_smiles_file)
        if 'Plant' in is_df.columns:
            for idx, row in is_df.iterrows():
                plant = row['Plant']
                ing = row['Ingredient name']
                if pd.notna(plant) and pd.notna(ing):
                    G.add_edge(plant, ing, weight=1.0, type='plant-ing')

except Exception as e:
    print(f"Error adding hierarchical edges: {e}")

# Calculate network metrics (only if graph not empty)
# NetworkX version compatibility: is_empty() might not exist in older/newer versions, use len(G.nodes()) > 0
if len(G.nodes()) > 0:
    degree_centrality = nx.degree_centrality(G)
    closeness_centrality = nx.closeness_centrality(G)
    betweenness_centrality = nx.betweenness_centrality(G)
    # eigenvector might fail on disconnected graph, use max_iter
    try:
        eigenvector_centrality = nx.eigenvector_centrality(G, max_iter=1000)
    except:
        eigenvector_centrality = {n: 0 for n in G.nodes()}

    # Combine all centrality metrics into a DataFrame
    centrality_measures = {
        'degree': degree_centrality,
        'closeness': closeness_centrality,
        'betweenness': betweenness_centrality,
        'eigenvector': eigenvector_centrality
    }

    df = pd.DataFrame(centrality_measures)

    # Calculate sum and sort by sum
    df['sum'] = df.sum(axis=1)
    df = df.sort_values('sum', ascending=False)

    # Save results as CSV file
    df.to_csv(f"{output_folder}/centrality_measures.csv")
    print(f"Network analysis complete. Nodes: {G.number_of_nodes()}, Edges: {G.number_of_edges()}")

    # Export for Cytoscape (Nodes and Edges)
    nx.write_graphml(G, f"{output_folder}/network_for_cytoscape.graphml")
    print(f"GraphML saved to {output_folder}/network_for_cytoscape.graphml")
else:
    print("Graph is empty. No network analysis performed.")


# --- Cell 51 ---
# Read CSV file
output_folder = "05_ppi"
csv_file = "05_ppi/centrality_measures.csv"
df = pd.read_csv(csv_file)

# Keep only Protein column
# df = df[["Protein"]]
df = df.iloc[:, 0:1]

# Keep only top 20 values
df = df.head(20)

os.makedirs(output_folder, exist_ok=True)

# Export as XLSX file
output_file = os.path.join(output_folder, "selected_proteins.xlsx")
df.to_excel(output_file, index=False)


# --- Cell 53: Skipped (Appears to be R code) ---
# --- Cell 54: Skipped (Appears to be R code) ---
# --- Cell 56 ---
# Skipped (R code)

# --- Cell 58 ---
# Skipped (R code)

# --- Cell 60 ---
# Skipped (R code)

# --- Cell 62 ---
# Skipped (R code)

# --- Cell 64 ---
# Skipped (R code)

# --- Cell 66 ---
# Skipped (R code)

# --- Cell 68 ---
# Skipped (R code)

# --- Cell 70 ---
# Skipped (R code)

# --- Cell 72 ---
# Skipped (R code)


# --- Cell 74: Skipped (Appears to be R code) ---
# --- Cell 75: Skipped (Appears to be R code) ---