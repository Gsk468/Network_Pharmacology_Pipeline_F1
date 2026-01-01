# Network Pharmacology Pipeline

A comprehensive workflow for network pharmacology analysis, integrating data processing, target prediction, disease gene association, protein-protein interaction (PPI) network construction, and enrichment analysis. The pipeline combines **Python** for data retrieval and processing with **R** for statistical analysis and visualization.

## Overview

This repository provides two ways to run the workflow:
1.  **Standalone Python Script (`Network_Pharmacology_Pipeline.py`)**: A fully automated pipeline suitable for batch processing and integration.
2.  **Jupyter Notebook (`Network_Pharmacology_Pipeline.ipynb`)**: An interactive environment for step-by-step execution and visualization.

## Key Features

*   **Automated Data Gathering**:
    *   **Drug Target Prediction**: Simulates target prediction (robust fallback) or integrates with prediction tools.
    *   **Disease Targets**: Fetches cancer-associated genes from the **OpenTargets Platform** API (default: Cancer, MONDO_0004992).
    *   **PPI Network**: Retrieves protein interactions from **STRING DB**.
*   **Analysis**:
    *   **Network Toplogy**: Calculates centrality measures (Degree, Closeness, Betweenness) to identify Hub Genes.
    *   **Enrichment**: Performs GO and KEGG pathway enrichment analysis using R (`clusterProfiler`).
*   **Visualization & Export**:
    *   Generates CSV files for **Cytoscape** visualization (`Drug-Target`, `PPI`, `Target-Pathway` networks).
    *   Produces static plots (Bar, Dot, Cnet, Heatmap) via R.

## Prerequisites

*   **Python 3.x**: Required for the main pipeline logic.
    *   Dependencies: `pandas`, `requests`, `networkx`, `matplotlib`, `pubchempy`, `chardet`, `openpyxl`.
*   **R**: Required for enrichment analysis.
    *   The pipeline attempts to automatically install required Bioconductor packages (`clusterProfiler`, `org.Hs.eg.db`, etc.) via `BiocManager`.

## Usage

### 1. Quick Start (Python Script)

The easiest way to run the analysis is using the standalone script. It handles input generation, API calls, and R execution automatically.

```bash
python Network_Pharmacology_Pipeline.py
```

**What happens:**
1.  **Input Generation**: If no input is found, it creates a default dataset of **10 known anticancer drugs** (e.g., Imatinib, Doxorubicin) in `01.Drug_Ingredients/05.merge/SMILES_qed0.67.csv`.
2.  **Pipeline Execution**:
    *   Predicts/Simulates targets for the input drugs.
    *   Fetches disease targets from OpenTargets.
    *   Intersects drug and disease targets.
    *   Constructs a PPI network via STRING DB.
    *   Calculates network centrality.
3.  **R Analysis**: Automatically invokes the R script to perform GO/KEGG enrichment.
4.  **Outputs**:
    *   **Cytoscape Files**: `07_figure/cytoscape_edges.csv`, `07_figure/cytoscape_nodes.csv`.
    *   **Enrichment Results**: `06_enrichment/`.
    *   **PPI Data**: `05_ppi/`.

### 2. Jupyter Notebook

For an interactive experience:

1.  Open `Network_Pharmacology_Pipeline.ipynb` in Jupyter (VSCode, JupyterLab, etc.).
2.  Run cells sequentially.
3.  The notebook contains the same robust logic as the script, allowing you to inspect intermediate dataframes and plots.

**Note on GitHub Viewing:**
GitHub often has issues rendering complex `.ipynb` files. You can view it using [nbviewer](https://nbviewer.org/) by pasting the repository URL.

## Directory Structure

*   `01.Drug_Ingredients/`: Input data (SMILES) and target prediction scripts.
*   `02_ingredients_targets/`: Intermediate target prediction results.
*   `03_diseases_targets/`: Disease targets fetched from OpenTargets.
*   `04_intersection_targets/`: Intersection of drug and disease targets.
*   `05_ppi/`: PPI network data, centrality measures, and selected hub genes.
*   `06_enrichment/`: GO and KEGG enrichment result tables.
*   `07_figure/`: Generated files for Cytoscape (Edges, Nodes) and R plots.
