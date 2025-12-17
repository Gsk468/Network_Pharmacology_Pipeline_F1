# Network_Pharmacology_Pipeline

A comprehensive workflow for network pharmacology analysis, focusing on the identification of active ingredients, target prediction, disease target interaction, and network analysis.

The pipeline involves:
1.  **Data Processing (Python):** Ingredient loading, target prediction (SwissTargetPrediction, SEA, PPB3), disease target acquisition (OpenTargets), and network construction.
2.  **Visualization & Enrichment (R):** GO/KEGG enrichment analysis and visualization.

**Note:** This workflow requires switching between Python and R kernels in Jupyter/VSCode or setting up a mixed environment.

## Environment Setup Guide

To ensure all scripts and the notebook run correctly, please follow these steps to set up your environment. Using **Conda** (Anaconda or Miniconda) is highly recommended to manage both Python and R dependencies.

### 1. Prerequisites

*   **Conda:** Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or Anaconda.
*   **Google Chrome:** Required for the target prediction script (`swiss-sea-ppb3-target-pred.py`) which uses Selenium. Ensure Google Chrome is installed on your system.

### 2. Create and Activate Conda Environment

Open your terminal or command prompt and run:

```bash
# Create a new environment named 'netpharm' with Python 3.9
conda create -n netpharm python=3.9 -y

# Activate the environment
conda activate netpharm
```

### 3. Install Python Dependencies

Install the required Python packages using pip:

```bash
pip install jupyter pandas numpy openpyxl requests tqdm beautifulsoup4 selenium webdriver-manager networkx matplotlib pubchempy chardet matplotlib-venn bioservices
```

*Note: `bioservices` might require additional build tools on some systems. If it fails, the pipeline has fallback logic, but full functionality is recommended.*

### 4. Install R and R Packages

The pipeline uses R for enrichment analysis. You can install R within the same Conda environment or use a system installation. Here we install via Conda for convenience:

```bash
# Install R base
conda install -c conda-forge r-base=4.2 -y

# Install R kernel for Jupyter
conda install -c conda-forge r-irkernel -y
```

**Installing R Libraries:**
Open an R terminal (type `R` in your console) or run these commands in an R script/notebook cell:

```R
# Install CRAN packages
install.packages(c("openxlsx", "ggplot2", "stringr", "GOplot", "ggnewscale", "circlize", "repr"), repos="http://cran.us.r-project.org")

# Install Bioconductor Manager
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos="http://cran.us.r-project.org")

# Install Bioconductor packages
BiocManager::install(c("enrichplot", "clusterProfiler", "DOSE", "topGO", "ComplexHeatmap", "org.Hs.eg.db"))
```

### 5. Running the Pipeline

1.  **Input Data:** Ensure `Sample-Plant-smiles.xlsx` is present in the root directory. This file should contain columns: `Phytochemical`, `Plant`, `SMILES`, `CID`.
2.  **Launch Jupyter:**
    ```bash
    jupyter notebook
    ```
3.  **Open Notebook:** Open `Network_Pharmacology_Pipeline.ipynb`.
4.  **Execute:** Run the cells sequentially.
    *   **Python Part:** Ensure the kernel is set to **Python 3 (netpharm)**.
    *   **R Part:** When reaching the Enrichment Analysis section (Step 6), switch the kernel to **R** (or ensure your environment supports R magics if configured). The notebook currently assumes manual kernel switching.

### 6. Troubleshooting

*   **Chrome/Selenium Issues:** If target prediction fails, ensure Google Chrome is installed. The script will try to use `webdriver-manager` to automatically download the matching ChromeDriver. If `headless` mode causes issues, you can modify the script arguments.
*   **Missing Dependencies:** If a package is missing, install it via `pip install <package_name>` (Python) or `install.packages("<package_name>")` (R).
