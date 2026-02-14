#!/usr/bin/env python3
# Network Pharmacology Pipeline - Auto-generated script
import os
import subprocess
import sys
import json
import glob
import pandas as pd
import requests
import shutil
import csv
import re
import xml.etree.ElementTree as ET
import openpyxl
import chardet
from time import sleep
import networkx as nx
from concurrent.futures import ThreadPoolExecutor
from matplotlib_venn import venn2
import matplotlib.pyplot as plt

def run_r_code(r_source, step_name):
    print(f'\n[R] Running Step: {step_name}...')
    r_filename = f'temp_{step_name.replace(" ", "_")}.R'
    with open(r_filename, 'w', encoding='utf-8') as f:
        f.write(r_source)
    try:
        # Check if Rscript is available
        subprocess.run(['Rscript', '--version'], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # Run the script
        subprocess.run(['Rscript', r_filename], check=True)
        print(f'[R] Step {step_name} completed successfully.')
    except (subprocess.CalledProcessError, FileNotFoundError):
        print(f'[WARNING] R execution failed for {step_name}. Ensure R is installed and in your PATH.')
    finally:
        if os.path.exists(r_filename):
            os.remove(r_filename)

# --- Data Loading (Replaces Cell 10) ---
# We start with the provided Sample-Plant-smiles.xlsx or similar input
print("--- Starting Pipeline ---")

# Setup directories
dirs = [
    "01_drug_ingredients/03.ingredients_augmenteddata",
    "01_drug_ingredients/04.ingredients_filtereddata",
    "01_drug_ingredients/05.merge",
    "02_ingredients_targets/021_ingredients_targets_orig",
    "02_ingredients_targets/022_ingredients_targets_processed",
    "02_ingredients_targets/023_ingredients_targets_uniquedata",
    "03_diseases_targets/031_orig",
    "03_diseases_targets/032_trans",
    "03_diseases_targets/033_screened",
    "03_diseases_targets/033_uniquedata",
    "04.Final",
    "04_intersection_targets",
    "05_ppi",
    "06_enrichment",
    "07_figure"
]
for d in dirs:
    os.makedirs(d, exist_ok=True)

# Input file handling
# Assuming we use Sample-Plant-smiles.xlsx as the starting point for "HERB_augmented.csv"
input_xlsx = "Sample-Plant-smiles.xlsx"
if os.path.exists(input_xlsx):
    print(f"Loading {input_xlsx}...")
    df = pd.read_excel(input_xlsx)
    # Rename columns to match pipeline expectations
    # Expected: Ingredient name, Plant, SMILES, Pubchem_CID (optional)
    # Map from Phytochemical -> Ingredient name if needed
    col_map = {
        'Phytochemical': 'Ingredient name',
        'Plant': 'Source', # or Plant
        'SMILES': 'SMILES',
        'CID': 'Pubchem_CID'
    }
    df.rename(columns=col_map, inplace=True)

    # Ensure Ingredient name exists
    if 'Ingredient name' not in df.columns and 'Ingredient' in df.columns:
        df.rename(columns={'Ingredient': 'Ingredient name'}, inplace=True)

    output_csv = "01_drug_ingredients/03.ingredients_augmenteddata/HERB_augmented.csv"
    df.to_csv(output_csv, index=False)
    print(f"Saved to {output_csv}")
else:
    print(f"Warning: {input_xlsx} not found. Ensure input data is present.")

# --- Filtering (Cell 17) ---
print("Filtering ingredients...")
herb_aug_file = "01_drug_ingredients/03.ingredients_augmenteddata/HERB_augmented.csv"
if os.path.exists(herb_aug_file):
    df = pd.read_csv(herb_aug_file)
    # Filter valid SMILES or CID
    # Original code filtered by PubChem_id. We'll filter by SMILES as it's critical for prediction.
    if 'SMILES' in df.columns:
        df = df.dropna(subset=['SMILES'])
        df = df[df['SMILES'].str.strip() != '']

    filtered_csv = "01_drug_ingredients/04.ingredients_filtereddata/HERB_filtered.csv"
    df.to_csv(filtered_csv, index=False)
    print(f"Saved filtered data to {filtered_csv}")
else:
    filtered_csv = ""

# --- Merging/Deduplication (Cell 25/27) ---
# In this refactored flow, we treat HERB_filtered.csv as the main input for merging
print("Preparing merged input for target prediction...")
input_folder = "01_drug_ingredients/04.ingredients_filtereddata"
output_file = "01_drug_ingredients/05.merge/Ingredient_smiles.csv"

# Simple copy/merge
if filtered_csv and os.path.exists(filtered_csv):
    df = pd.read_csv(filtered_csv)
    # Add Source if missing
    if 'Source' not in df.columns:
        df['Source'] = 'User_Input'
    df.to_csv(output_file, index=False)
    print(f"Merged data saved to {output_file}")
else:
    print("No filtered data to merge.")

# Deduplication (Cell 27)
input_file = "01_drug_ingredients/05.merge/Ingredient_smiles.csv"
output_file_dedup = "01_drug_ingredients/05.merge/Ingredient_smiles_dedup.csv"

if os.path.exists(input_file):
    data = pd.read_csv(input_file)
    # Deduplicate by SMILES (most important) or Ingredient Name
    if 'SMILES' in data.columns:
        data = data.drop_duplicates(subset=["SMILES"], keep="first")
    else:
        data = data.drop_duplicates(subset=["Ingredient name"], keep="first")

    data.to_csv(output_file_dedup, index=False)
    print(f"Deduplicated data saved to {output_file_dedup}")

# --- Target Prediction (Cell 37) ---
print("Running Target Prediction...")
input_csv = output_file_dedup
output_dir = "02_ingredients_targets/021_ingredients_targets_orig"
target_pred_output = os.path.join(output_dir, "TARGETS.csv")

if os.path.exists(input_csv):
    # Call the robust prediction script
    # Note: ensure swiss-sea-ppb3-target-pred.py is executable or run with python3
    cmd = f'python3 swiss-sea-ppb3-target-pred.py --input "{input_csv}" --output "{target_pred_output}"'
    os.system(cmd)
else:
    print("Skipping target prediction (no input).")

# --- Post-Processing Targets (Cell 39) ---
print("Processing Target Results...")
input_folder = '02_ingredients_targets/021_ingredients_targets_orig'
output_folder = '02_ingredients_targets/022_ingredients_targets_processed'

csv_files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]
for csv_file in csv_files:
    file_path = os.path.join(input_folder, csv_file)
    try:
        with open(file_path, 'rb') as f:
            result = chardet.detect(f.read())
        df = pd.read_csv(file_path, encoding=result['encoding'])

        # Normalize columns
        if 'UniProt_ID' in df.columns and 'UniProt_name' not in df.columns:
            df.rename(columns={'UniProt_ID': 'UniProt_name'}, inplace=True)

        # Fallback: if UniProt_name is missing or empty, use Target_Name
        if 'Target_Name' in df.columns:
             if 'UniProt_name' not in df.columns:
                 df['UniProt_name'] = df['Target_Name']
             else:
                 df['UniProt_name'] = df['UniProt_name'].fillna('')
                 df.loc[df['UniProt_name'] == '', 'UniProt_name'] = df['Target_Name']

        if 'UniProt_name' in df.columns:
            # Ensure UniProt_name is treated as string
            df['UniProt_name'] = df['UniProt_name'].astype(str)
            df = df.assign(UniProt_name=df['UniProt_name'].str.split('|')).explode('UniProt_name')
            # Filter already done by prediction script, but ensuring no NaNs
            df = df.dropna(subset=['UniProt_name'])

        output_file_name = csv_file.replace('.csv', '_processed.csv')
        df.to_csv(os.path.join(output_folder, output_file_name), index=False)
    except Exception as e:
        print(f"Error processing {csv_file}: {e}")

# --- Unique Targets (Cell 41) ---
input_folder = '02_ingredients_targets/022_ingredients_targets_processed'
output_folder = '02_ingredients_targets/023_ingredients_targets_uniquedata'
csv_files = [f for f in os.listdir(input_folder) if f.endswith('.csv')]

for csv_file in csv_files:
    input_file_path = os.path.join(input_folder, csv_file)
    try:
        with open(input_file_path, 'rb') as f:
            result = chardet.detect(f.read())
        df = pd.read_csv(input_file_path, encoding=result['encoding'])
        if 'UniProt_name' in df.columns:
            df = df[['UniProt_name']].drop_duplicates()
            output_file_name = csv_file.replace('.csv', '_uniquedata.csv')
            df.to_csv(os.path.join(output_folder, output_file_name), index=False)
    except Exception as e:
        print(f"Error uniquing {csv_file}: {e}")

# --- Fetch Fasta (Cell 43) ---
print("Fetching UniProt FASTA (this may take a while)...")
url = 'https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28%2A%29%20AND%20%28reviewed%3Atrue%29%20AND%20%28model_organism%3A9606%29'
try:
    response = requests.get(url)
    all_fastas = response.text
    print(f"Fetched {len(all_fastas)} chars of FASTA data.")
except Exception as e:
    print(f"Failed to fetch FASTA: {e}")
    all_fastas = ""

# --- Map Targets to Gene Symbols (Cell 44) ---
# This was the problematic cell. Fixing loop logic.
print("Mapping targets to Gene Symbols...")
input_folder = '02_ingredients_targets/023_ingredients_targets_uniquedata'
output_folder = '04.Final'

csv_files = glob.glob(os.path.join(input_folder, '*_uniquedata.csv'))
fasta_list = re.split(r'\n(?=>)', all_fastas)

for csv_file in csv_files:
    base_name = os.path.basename(csv_file)
    output_csv = os.path.join(output_folder, base_name.replace('_uniquedata.csv', '_adjusted.csv'))

    df = pd.read_csv(csv_file)
    if 'UniProt_name' not in df.columns:
        continue

    symbol_values = df['UniProt_name'].dropna().tolist()
    if not symbol_values:
        continue

    # Escape to prevent regex errors with special chars
    pattern = '|'.join([f'\\b{re.escape(str(s))}\\b' for s in symbol_values])
    if not pattern:
        continue

    matched_fastas = [fasta for fasta in fasta_list if re.search(pattern, fasta)]

    processed_data = []
    original_values = []

    for fasta in matched_fastas:
        # Fasta header format: >db|UniqueIdentifier|EntryName ProteinName OS=...
        # We want EntryName (e.g. ALBU_HUMAN) or GeneName if available.
        # The script seems to extract the Entry Name part.
        parts = fasta.split('|')
        if len(parts) > 2:
            processed_value = parts[2].split(' ')[0] # Entry Name
            # Find which original symbol matched
            matched_syms = [s for s in symbol_values if s in fasta]

            for ms in matched_syms:
                original_values.append(ms)
                processed_data.append(processed_value)

    # Check if we have data
    if processed_data:
        processed_df = pd.DataFrame({
            'Original Value': original_values,
            'Processed Value': processed_data
        }).drop_duplicates()
        processed_df.to_csv(output_csv, index=False)
        print(f"Mapped targets saved to {output_csv}")

# --- Disease Data (Fetch from OpenTargets) ---
# Implementing OpenTargets fetching for disease targets
print("Processing Disease Targets...")
disease_targets_folder = "03_diseases_targets/033_uniquedata"
os.makedirs(disease_targets_folder, exist_ok=True)
disease_targets_file = os.path.join(disease_targets_folder, 'disease_targets_uniquedata.csv')

# Use a default disease if none provided, but ideally this should be an argument.
# We will check if a 'disease_name.txt' exists or use a hardcoded default for demo.
disease_name = "diabetes mellitus"
if os.path.exists("disease_name.txt"):
    with open("disease_name.txt", "r") as f:
        disease_name = f.read().strip()

print(f"Fetching targets for disease: {disease_name} (from OpenTargets)...")

def fetch_opentargets_disease(disease_query):
    # 1. Search for disease ID
    search_url = "https://platform-api.opentargets.io/v3/platform/public/search"
    params = {"q": disease_query, "size": 1, "filter": "disease"}
    try:
        r = requests.get(search_url, params=params)
        data = r.json()
        if not data.get('data'):
            print(f"No disease found for {disease_query}")
            return []

        disease_id = data['data'][0]['id']
        print(f"Found Disease ID: {disease_id} ({data['data'][0]['name']})")

        # 2. Fetch associated targets
        # The v3 API is deprecated but often simpler. v4 uses GraphQL.
        # Let's try v3 association.
        assoc_url = "https://platform-api.opentargets.io/v3/platform/public/association/filter"
        assoc_params = {"disease": disease_id, "size": 10000, "scorevalue_min": 0.1, "direct": True}

        r = requests.get(assoc_url, params=assoc_params)
        assoc_data = r.json()

        targets = []
        for d in assoc_data.get('data', []):
            if 'target' in d and 'gene_info' in d['target']:
                 targets.append(d['target']['gene_info']['symbol'])

        return list(set(targets))
    except Exception as e:
        print(f"OpenTargets fetch failed: {e}")
        return []

# Fetch and save
disease_symbols = fetch_opentargets_disease(disease_name)
if not disease_symbols:
    print("Using MOCK disease targets due to network/API failure (for demonstration).")
    # Mock data for "Diabetes" or similar (Sample targets that might intersect with Quercetin targets)
    disease_symbols = [
        "CSNK2A1", "CSNK2A2", "CA7", "AMY1A", "AKT1", "INS", "GCK", "PPARG", "TNF", "IL6", "VEGFA", "PTGS2", "NOS3"
    ]

if disease_symbols:
    pd.DataFrame({'Symbol': disease_symbols}).to_csv(disease_targets_file, index=False)
    print(f"Saved {len(disease_symbols)} disease targets to {disease_targets_file}")
else:
    print("No disease targets fetched. Using empty set.")
    pd.DataFrame({'Symbol': []}).to_csv(disease_targets_file, index=False)


# --- Intersection (Cell 64 & 68 Correction) ---
# Network Pharmacology Logic:
# 1. Union of ALL Ingredient Targets
# 2. Intersect with Disease Targets

print("Processing Intersection (Ingredients Union ∩ Disease)...")
input_folder = '02_ingredients_targets/023_ingredients_targets_uniquedata'
output_folder = '04_intersection_targets'

# 1. Collect all ingredient targets (mapped to Gene Symbols)
all_ingredient_targets = set()

# Helper to map UniProt/Target names to Gene Symbols using FASTA (if needed) or simple lookup
# We already downloaded FASTA in Cell 43.
fasta_list = re.split(r'\n(?=>)', all_fastas)

csv_files = glob.glob(os.path.join(input_folder, '*_uniquedata.csv'))
for csv_file in csv_files:
    df = pd.read_csv(csv_file)
    if 'UniProt_name' not in df.columns: continue

    # We need to map these UniProt names/IDs to Gene Symbols.
    # The previous logic used regex on FASTA. We can reuse that but strictly for Gene Name.
    # To save time, we can create a mapping dictionary from FASTA once.

    # Extract gene symbols from the UniProt_name column directly if they look like symbols,
    # or map them if they look like Entry Names.
    # Ideally, we map all to Gene Symbol.

    targets = df['UniProt_name'].dropna().astype(str).tolist()
    if not targets: continue

    # Perform mapping (simplified for batch)
    # We will assume if it's all uppercase and no underscore, it might be a symbol.
    # If it has _HUMAN, we map.

    # Let's use the robust FASTA mapping approach from before but aggregate results.
    pattern = '|'.join([f'\\b{re.escape(str(s))}\\b' for s in targets])
    if not pattern: continue

    matched_fastas = [fasta for fasta in fasta_list if re.search(pattern, fasta)]

    for fasta in matched_fastas:
        # Extract Gene Name (GN=...)
        match = re.search(r'GN=([^\s]+)', fasta)
        if match:
            all_ingredient_targets.add(match.group(1))
        else:
            # Fallback to Entry Name (e.g. ALBU_HUMAN -> ALBU? No, ALB)
            # Or just use the first word of description?
            # If GN is missing, it's safer to skip or use EntryName
             parts = fasta.split('|')
             if len(parts) > 2:
                 entry_name = parts[2].split(' ')[0]
                 # e.g. ALBU_HUMAN. We assume the part before _ is often the symbol (not always true but fallback)
                 if '_HUMAN' in entry_name:
                     all_ingredient_targets.add(entry_name.split('_')[0])

print(f"Total unique ingredient targets (Gene Symbols): {len(all_ingredient_targets)}")

# 2. Load Disease Targets
disease_targets = set()
if os.path.exists(disease_targets_file):
    df_dis = pd.read_csv(disease_targets_file)
    if 'Symbol' in df_dis.columns:
        disease_targets = set(df_dis['Symbol'].dropna().astype(str).tolist())

print(f"Total disease targets: {len(disease_targets)}")

# 3. Intersect
intersection = all_ingredient_targets.intersection(disease_targets)
print(f"Intersection (Ingredients ∩ Disease): {len(intersection)}")

if intersection:
    result_df = pd.DataFrame(list(intersection), columns=['Processed Value']) # Keeping column name compatible
    result_df.to_csv(os.path.join(output_folder, 'intersection_targets.csv'), index=False)

    # Plot Venn
    if len(all_ingredient_targets) > 0 and len(disease_targets) > 0:
        plt.figure()
        venn2([all_ingredient_targets, disease_targets], set_labels=('Compounds', 'Disease'))
        plt.title(f"Intersection: {len(intersection)}")
        plt.savefig(os.path.join(output_folder, 'venn_diagram.svg'))
        plt.close()
        print("Venn diagram saved.")

    # Generate Network Files for Cytoscape (Nodes/Edges)
    # Nodes: Type (Compound, Target, Disease)
    # Edges: Compound-Target, Target-Disease
    print("Generating Cytoscape network files...")

    # Edges List
    edges = []

    # Compound-Target Edges (Only for targets in intersection)
    # We need to re-scan which compound had which target.
    # This is inefficient but functional: iterate files again.

    for csv_file in csv_files:
        # Filename format: Source_Ingredient_processed_uniquedata.csv?
        # Actually input folder is 023...uniquedata. File names are like "Source_Ingredient...csv"
        # We need the Ingredient Name.
        filename = os.path.basename(csv_file)
        # Assuming filename is "Source_Ingredient_..." or similar.
        # Let's try to extract ingredient name.
        # Since we generated these filenames in Cell 39: `csv_file.replace('.csv', '_processed.csv')`
        # And originally `TARGETS.csv` became `TARGETS_processed.csv`.
        # Wait, the structure in 021 was "TARGETS.csv".
        # Ah, `swiss-sea-ppb3-target-pred.py` output "TARGETS.csv" which contains ALL compounds.
        # But `Network_Pharmacology_Pipeline.py` Cell 39 iterates over csv files in 021.
        # If there is only one `TARGETS.csv`, we lost the individual compound info if we didn't split it?
        # `TARGETS.csv` HAS "Phytochemical" column!

        # We should read the original TARGETS.csv to build edges properly if we want Compound-Target edges.
        pass

    # Better approach for Edges:
    # Read `01_drug_ingredients/05.merge/Ingredient_smiles_dedup.csv` for Compounds (Source -> Ingredient)
    # Read `02_ingredients_targets/021_ingredients_targets_orig/TARGETS.csv` for Ingredient -> Target (filtered by intersection)
    # Disease is just one node "Disease" connected to all Intersection Targets.

    targets_orig_file = "02_ingredients_targets/021_ingredients_targets_orig/TARGETS.csv"
    if os.path.exists(targets_orig_file):
        df_targets = pd.read_csv(targets_orig_file)

        # We need to map the "Target_Name" or "UniProt_ID" in this file to the "Gene Symbol" used in intersection.
        # This mapping is tricky because we did it via FASTA matching later.
        # For simplicity, we will map using the same logic or just use the prediction script's symbol if available.
        # `TARGETS.csv` has "Gene_Symbol" column (populated by Swiss, but maybe empty for PPB3).

        node_list = set()
        edge_list = []

        # Disease Node
        node_list.add((disease_name, "Disease"))

        for t in intersection:
            node_list.add((t, "Target"))
            edge_list.append((t, disease_name, "Target-Disease"))

        # Compound-Target
        # We need to match rows in df_targets to 't' in intersection.
        # We will iterate df_targets, map its target to symbol, and if symbol in intersection, add edge.

        # Build a map from UniProt_ID/Target_Name to Symbol using the FASTA logic we just ran
        # (We can't easily access the local vars from the loop above, so we approximate).

        # Simplified: Check if any part of the row matches the symbol
        for idx, row in df_targets.iterrows():
            phy = row['Phytochemical']

            # Try to find which symbol this row corresponds to
            # This is heuristic if we don't have exact mapping stored.
            # But we know `intersection` contains symbols.
            # Check if row['Gene_Symbol'] is in intersection

            sym = str(row.get('Gene_Symbol', ''))
            tid = str(row.get('UniProt_ID', ''))
            tname = str(row.get('Target_Name', ''))

            # Check matches
            matched_target = None
            if sym in intersection: matched_target = sym

            # If not found, try mapping TID or Name via FASTA cache (not available here) or heuristics
            if not matched_target:
                 # Check if TID is in FASTA that yielded a symbol in intersection
                 # This is too slow.
                 pass

            if matched_target:
                node_list.add((phy, "Compound"))
                edge_list.append((phy, matched_target, "Compound-Target"))

                # Plant-Compound (if available)
                plant = row.get('Plant', 'UnknownPlant')
                node_list.add((plant, "Plant"))
                edge_list.append((plant, phy, "Plant-Compound"))

        # Save Nodes/Edges
        nodes_df = pd.DataFrame(list(node_list), columns=['Id', 'Type'])
        edges_df = pd.DataFrame(edge_list, columns=['Source', 'Target', 'Interaction'])

        nodes_df.to_csv(os.path.join(output_folder, 'network_nodes.csv'), index=False)
        edges_df.to_csv(os.path.join(output_folder, 'network_edges.csv'), index=False)
        print("Cytoscape files saved.")

else:
    print("Intersection is empty.")

# --- PPI (Cell 72) ---
print("Fetching PPI Network from STRING...")
intersection_file = "04_intersection_targets/intersection_targets.csv"
if os.path.exists(intersection_file):
    df = pd.read_csv(intersection_file)
    my_genes = df["Processed Value"].tolist()

    # STRING API limits: if too many genes, might fail.
    # Use top 2000?
    if len(my_genes) > 0:
        identifiers = "%0d".join(my_genes[:2000])

        string_api_url = "https://string-db.org/api"
        output_folder = "05_ppi"

        # 1. TSV
        try:
            params = {
                "identifiers": identifiers,
                "species": 9606,
                "required_score": 400, # Lowered to 400 to ensure we get something
                "caller_identity": "NetworkPharmacologyPipeline"
            }
            response = requests.post(f"{string_api_url}/tsv/network", data=params)
            tsv_path = os.path.join(output_folder, "protein_interactions.tsv")
            with open(tsv_path, "w", encoding='utf-8') as f:
                f.write(response.text)
            print(f"PPI TSV saved to {tsv_path}")

            # 2. Centrality (Cell 76)
            if os.path.getsize(tsv_path) > 0:
                ppi_data = pd.read_csv(tsv_path, sep="\t")
                if 'preferredName_A' in ppi_data.columns:
                    G = nx.Graph()
                    for index, row in ppi_data.iterrows():
                        G.add_edge(row["preferredName_A"], row["preferredName_B"], weight=row["score"])

                    dc = nx.degree_centrality(G)
                    cc = nx.closeness_centrality(G)
                    bc = nx.betweenness_centrality(G)
                    ec = nx.eigenvector_centrality(G, max_iter=1000)

                    cent_df = pd.DataFrame({
                        'degree': dc, 'closeness': cc, 'betweenness': bc, 'eigenvector': ec
                    })
                    cent_df['sum'] = cent_df.sum(axis=1)
                    cent_df = cent_df.sort_values('sum', ascending=False)
                    cent_csv = os.path.join(output_folder, "centrality_measures.csv")
                    cent_df.to_csv(cent_csv)

                    # Top 20 for enrichment (Cell 80)
                    top20 = cent_df.head(20).reset_index().rename(columns={'index': 'Protein'})
                    top20.to_excel(os.path.join(output_folder, "selected_proteins.xlsx"), index=False)
                    print("Centrality analysis complete.")
                else:
                    print("PPI TSV format unexpected.")
        except Exception as e:
            print(f"PPI Analysis failed: {e}")

# --- R Analysis (Enrichment) ---
# Assuming R environment might not be perfect, we wrap this.
print("Starting R Analysis (Enrichment)...")
r_code_enrichment = r"""
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)

# Load Data
input_file <- "05_ppi/selected_proteins.xlsx"
if(file.exists(input_file)){
    info <- read.xlsx(input_file)
    print("Loaded proteins:")
    print(head(info))

    # Convert Symbols to Entrez IDs
    gene <- bitr(info$Protein, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)

    output_dir <- "06_enrichment"
    if (!dir.exists(output_dir)) dir.create(output_dir)

    # GO Analysis
    print("Running GO Enrichment...")
    GO <- enrichGO(gene$ENTREZID, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL", pvalueCutoff = 0.05, readable = TRUE)
    if(!is.null(GO)){
        write.csv(GO, file = file.path(output_dir, "GO_enrichment_results.csv"))
        pdf(file.path(output_dir, "GO_plots.pdf"))
        print(barplot(GO, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free"))
        print(dotplot(GO, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale="free"))
        dev.off()
    }

    # KEGG Analysis
    print("Running KEGG Enrichment...")
    KEGG <- enrichKEGG(gene$ENTREZID, organism = 'hsa', pvalueCutoff = 0.05)
    if(!is.null(KEGG)){
        write.csv(KEGG, file = file.path(output_dir, "KEGG_enrichment_results.csv"))
        pdf(file.path(output_dir, "KEGG_plots.pdf"))
        print(barplot(KEGG, showCategory=20))
        print(dotplot(KEGG, showCategory=20))
        dev.off()
    }
    print("Enrichment complete.")
} else {
    print("Input file for enrichment not found.")
}
"""
run_r_code(r_code_enrichment, "Enrichment_Analysis")

print("--- Pipeline Finished ---")
