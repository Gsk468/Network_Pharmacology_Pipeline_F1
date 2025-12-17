#!/usr/bin/env python3

"""

Combined Target Prediction Tool (FIXED VERSION - FORMATTING CORRECTED)

This script integrates SwissTargetPrediction and SEA (Similarity Ensemble Approach)

to perform comprehensive target prediction for a list of SMILES strings.

IMPORTANT FIX: Human target filtering now checks the Target_Key (UniProt ID) column

for "_HUMAN" suffix instead of relying on Description field, which ensures accurate

filtering of human targets.

FORMATTING FIX: Combined results now have proper column arrangement and consistent

ranking per SMILES per database.

Key features:

- Processes SMILES from an input file

- Runs predictions on both SwissTargetPrediction and SEA databases

- Filters for HUMAN targets only (by checking UniProt ID suffix)

- Saves separate CSV files for each database

- Saves a combined CSV file with all predictions (PROPERLY FORMATTED)

- Recalculates ranking per SMILES per database for consistency

- Configurable headless browser operation

- Robust error handling and rate limiting

"""

import time

import re

import random

import argparse

import pandas as pd

from pathlib import Path

from selenium import webdriver

from selenium.webdriver.chrome.service import Service

from selenium.webdriver.chrome.options import Options

from selenium.webdriver.common.by import By

from selenium.webdriver.support.ui import WebDriverWait

from selenium.webdriver.support import expected_conditions as EC

from selenium.common.exceptions import TimeoutException

from webdriver_manager.chrome import ChromeDriverManager

# Configuration

SWISSTARGET_URL = "https://www.swisstargetprediction.ch"

SEA_URL = "https://sea.bkslab.org/"

USER_AGENTS = [

"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/109.0.0.0 Safari/537.36",

"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/109.0.0.0 Safari/537.36",

"Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/108.0.0.0 Safari/537.36",

"Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/16.1 Safari/605.1.15",

"Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/108.0.0.0 Safari/537.36",

]


def setup_driver(headless=False):

    """Configures and returns a Selenium Chrome WebDriver."""

    opts = Options()

    user_agent = random.choice(USER_AGENTS)

    opts.add_argument(f"user-agent={user_agent}")

    print(f"Using User-Agent: {user_agent}")

    opts.add_argument("--no-sandbox")

    opts.add_argument("--disable-dev-shm-usage")

    opts.add_argument("--disable-gpu")

    opts.add_argument("--window-size=1366,768")

    opts.add_experimental_option("excludeSwitches", ["enable-automation"])

    opts.add_experimental_option("useAutomationExtension", False)

    if headless:

        opts.add_argument("--headless=new")

    service = Service(ChromeDriverManager().install())

    driver = webdriver.Chrome(service=service, options=opts)

    driver.set_page_load_timeout(300)

    driver.execute_script(

        "Object.defineProperty(navigator, 'webdriver', {get: () => undefined})"

    )

    return driver


def human_like_typing(element, text):

    """Simulates human-like typing into a web element."""

    for char in text:

        element.send_keys(char)

        time.sleep(random.uniform(0.05, 0.2))


# ==================== SWISSTARGETPREDICTION FUNCTIONS ====================

def swisstarget_submit(driver, smiles, species="Homo sapiens"):

    """Submits a SMILES string to SwissTargetPrediction."""

    print(f"[SwissTarget] Submitting: {smiles}")

    driver.get(SWISSTARGET_URL + "/")

    WebDriverWait(driver, 30).until(

        EC.presence_of_element_located((By.TAG_NAME, "body"))

    )

    time.sleep(random.uniform(1, 3))

    # Input SMILES

    smiles_input = WebDriverWait(driver, 20).until(

        EC.presence_of_element_located((By.CSS_SELECTOR, "input[name='smiles']"))

    )

    smiles_input.clear()

    human_like_typing(smiles_input, smiles)

    time.sleep(random.uniform(0.5, 1.5))

    # Select species (Human only)

    try:

        species_sel = driver.find_element(

            By.CSS_SELECTOR, "select[name='species'], select[name='organism']"

        )

        species_sel.send_keys(species)

        time.sleep(random.uniform(0.5, 1.5))

    except Exception:

        print("        [SwissTarget] Species selector not found, proceeding...")

    # Submit

    submit_button_selectors = ["input[type='submit']", "input[value*='Predict']", "#submitbutton"]

    submit = None

    for sel in submit_button_selectors:

        try:

            submit = driver.find_element(By.CSS_SELECTOR, sel)

            break

        except Exception:

            continue

    if submit:

        driver.execute_script("arguments[0].click();", submit)

    else:

        raise RuntimeError("[SwissTarget] Could not find submit button")


def swisstarget_get_results(driver, max_wait=300):

    """Waits for and returns the SwissTarget results table."""

    print("[SwissTarget] Waiting for results...")

    def find_table(drv):

        tables = drv.find_elements(By.TAG_NAME, "table")

        for t in tables:

            headers = [th.text.strip().lower() for th in t.find_elements(By.TAG_NAME, "th")]

            header_line = " ".join(headers)

            if ("target" in header_line and "uniprot id" in header_line and

                "probability" in header_line and "common" in header_line):

                return t

        return False

    return WebDriverWait(driver, max_wait).until(find_table)


def swisstarget_parse_table(table, smiles):

    """Parses SwissTarget results table."""

    rows = []

    headers = [th.text.strip().lower() for th in table.find_elements(By.TAG_NAME, "th")]

    header_map = {

        'target': next((h for h in headers if 'target' in h), None),

        'gene': next((h for h in headers if 'common' in h), None),

        'uniprot': next((h for h in headers if 'uniprot' in h), None),

        'chembl': next((h for h in headers if 'chembl' in h), None),

        'class': next((h for h in headers if 'class' in h), None),

        'prob': next((h for h in headers if 'probability' in h), None),

    }

    try:

        col_indices = {key: headers.index(val) for key, val in header_map.items() if val}

    except ValueError as e:

        raise RuntimeError(f"[SwissTarget] Missing required columns: {e}")

    prob_regex = re.compile(r"[0-9]*\.?[0-9]+")

    trs = table.find_elements(By.TAG_NAME, "tr")

    for i, tr in enumerate(trs[1:], start=1):

        tds = tr.find_elements(By.TAG_NAME, "td")

        if len(tds) < len(col_indices):

            continue

        def get_text(col_name):

            idx = col_indices.get(col_name)

            return tds[idx].text.strip() if idx is not None and idx < len(tds) else ""

        target_name = get_text('target')

        uniprot_id = get_text('uniprot')

        if not target_name or not uniprot_id:

            continue

        prob_text = get_text('prob')

        m = prob_regex.search(prob_text)

        prob_val = float(m.group(0)) if m else 0.0

        rows.append({

            "SMILES": smiles,

            "Database": "SwissTargetPrediction",

            "Rank": i,

            "Target_Name": target_name,

            "Gene_Symbol": get_text('gene'),

            "UniProt_ID": uniprot_id,

            "ChEMBL_ID": get_text('chembl'),

            "Target_Class": get_text('class'),

            "Target_Key": uniprot_id,

            "Probability": prob_val,

            "P_Value": "",

            "Max_Tc": "",

        })

    return rows


# ==================== SEA FUNCTIONS ====================

def sea_submit(driver, smiles):

    """Submits a SMILES string to SEA."""

    print(f"[SEA] Submitting: {smiles}")

    driver.get(SEA_URL + "/")

    WebDriverWait(driver, 30).until(

        EC.presence_of_element_located((By.TAG_NAME, "body"))

    )

    time.sleep(random.uniform(1, 3))

    # Input SMILES

    smiles_input = WebDriverWait(driver, 20).until(

        EC.presence_of_element_located((By.CSS_SELECTOR, "input[placeholder='Paste SMILES or try the example below']"))

    )

    smiles_input.clear()

    human_like_typing(smiles_input, smiles)

    time.sleep(random.uniform(0.5, 1.5))

    # Submit

    submit_button = WebDriverWait(driver, 10).until(

        EC.element_to_be_clickable((By.XPATH, "//button[contains(., 'Try SEA')]"))

    )

    submit_button.click()


def sea_get_results(driver, max_wait=300):

    """Waits for and returns the SEA results table."""

    print("[SEA] Waiting for results...")

    WebDriverWait(driver, max_wait).until(EC.url_contains("/jobs/"))

    table = WebDriverWait(driver, 60).until(

        EC.presence_of_element_located((By.CSS_SELECTOR, "table.table.table-bordered"))

    )

    return table


def is_human_target(target_key):

    """

    FIXED: Checks if a target is human by examining the UniProt ID (Target_Key).

    Human targets end with "_HUMAN" suffix.

    Examples:

    - CP2CJ_HUMAN -> Human

    - ANDR_HUMAN -> Human

    - KCNJ3_HUMAN -> Human

    - RORC_MOUSE -> Not human

    - MAOB_MOUSE -> Not human

    - CHIB1_ASPFM -> Not human

    """

    if not target_key or not isinstance(target_key, str):

        return False

    # Check if Target_Key ends with _HUMAN

    return target_key.upper().endswith("_HUMAN")


def sea_parse_table(table, smiles):

    """

    Parses SEA results table and filters for HUMAN targets only.

    FIXED: Now checks Target_Key (UniProt ID) column for "_HUMAN" suffix

    instead of relying on Description field which may be incomplete or inconsistent.

    """

    rows = []

    thead_ths = table.find_elements(By.TAG_NAME, "thead")[0].find_elements(By.TAG_NAME, "th")

    thead_header_texts = [th.text.strip().lower() for th in thead_ths]

    col_mapping = {

        "Target_Key": "target key",

        "Target_Name": "target name",

        "Description": "description",

        "P_Value": "p-value",

        "Max_Tc": "maxtc",

    }

    actual_td_indices = {}

    for output_col, header_text in col_mapping.items():

        try:

            actual_td_indices[output_col] = thead_header_texts.index(header_text) - 1

        except ValueError:

            raise RuntimeError(f"[SEA] Missing expected header: '{header_text}'")

    tbody_rows = table.find_elements(By.TAG_NAME, "tbody")[0].find_elements(By.TAG_NAME, "tr")

    rank = 0

    for tr in tbody_rows:

        tr_class = tr.get_attribute("class")

        if tr_class and "spanning info" in tr_class:

            continue

        rank += 1

        tds = tr.find_elements(By.TAG_NAME, "td")

        if len(tds) < len(col_mapping):

            continue

        target_key = tds[actual_td_indices["Target_Key"]].text.strip()

        target_name = tds[actual_td_indices["Target_Name"]].text.strip()

        description = tds[actual_td_indices["Description"]].text.strip()

        # FIXED: Check UniProt ID (Target_Key) for _HUMAN suffix

        if not is_human_target(target_key):

            print(f"[SEA] Skipping non-human target: {target_name} | {target_key} | {description}")

            continue

        rows.append({

            "SMILES": smiles,

            "Database": "SEA",

            "Rank": rank,

            "Target_Name": target_name,

            "Gene_Symbol": "",

            "UniProt_ID": target_key,

            "ChEMBL_ID": "",

            "Target_Class": description,

            "Target_Key": target_key,

            "Probability": "",

            "P_Value": tds[actual_td_indices["P_Value"]].text.strip(),

            "Max_Tc": tds[actual_td_indices["Max_Tc"]].text.strip(),

        })

    return rows


# ==================== FORMATTING FUNCTIONS ====================

def format_combined_results(swisstarget_rows, sea_rows, output_file):

    """

    Formats combined results with proper column arrangement and ranking.

    FIXED: Recalculates rank per SMILES per database and ensures

    consistent column arrangement.

    """

    all_rows = swisstarget_rows + sea_rows

    # Create DataFrame

    df_combined = pd.DataFrame(all_rows)

    # Define proper column order

    columns_order = [

        "SMILES",

        "Database",

        "Target_Name",

        "Gene_Symbol",

        "UniProt_ID",

        "ChEMBL_ID",

        "Target_Class",

        "Target_Key",

        "Probability",

        "P_Value",

        "Max_Tc",

        "Rank"

    ]

    # Reorder columns

    df_combined = df_combined[columns_order]

    # Recalculate rank per SMILES per Database

    df_combined['Rank'] = df_combined.groupby(['SMILES', 'Database']).cumcount() + 1

    # Sort by SMILES, Database, then Rank

    df_combined = df_combined.sort_values(['SMILES', 'Database', 'Rank']).reset_index(drop=True)

    # Fill NaN values with empty strings for cleaner CSV output

    df_combined = df_combined.fillna('')

    # Save to CSV

    df_combined.to_csv(output_file, index=False)

    print(f"\n[Combined] Saved {len(df_combined)} rows to {output_file}")

    print(f"[Combined] Columns: {', '.join(columns_order)}")

    return df_combined


# ==================== MAIN WORKFLOW ====================

def run_combined_prediction(smiles_list, output_dir="results", headless=False):

    """Runs target prediction on both databases."""

    output_dir = Path(output_dir)

    output_dir.mkdir(exist_ok=True)

    swisstarget_rows = []

    sea_rows = []

    driver = setup_driver(headless=headless)

    try:

        for idx, smi in enumerate(smiles_list, 1):

            print(f"\n{'='*60}")

            print(f"Processing [{idx}/{len(smiles_list)}]: {smi}")

            print(f"{'='*60}")

            # SwissTargetPrediction

            print("\n--- Running SwissTargetPrediction ---")

            try:

                swisstarget_submit(driver, smi, species="Homo sapiens")

                table = swisstarget_get_results(driver)

                rows = swisstarget_parse_table(table, smi)

                print(f"[SwissTarget] Successfully parsed {len(rows)} human targets")

                swisstarget_rows.extend(rows)

            except Exception as e:

                print(f"[SwissTarget] Error for {smi}: {e}")

                swisstarget_rows.append({

                    "SMILES": smi, "Database": "SwissTargetPrediction", "Rank": 0,

                    "Target_Name": "Processing failed", "Gene_Symbol": "",

                    "UniProt_ID": "", "ChEMBL_ID": "", "Target_Class": "",

                    "Target_Key": "", "Probability": 0.0, "P_Value": "", "Max_Tc": "",

                })

            time.sleep(random.uniform(3, 6))

            # SEA

            print("\n--- Running SEA ---")

            try:

                sea_submit(driver, smi)

                table = sea_get_results(driver)

                rows = sea_parse_table(table, smi)

                print(f"[SEA] Successfully parsed {len(rows)} human targets")

                sea_rows.extend(rows)

            except Exception as e:

                print(f"[SEA] Error for {smi}: {e}")

                sea_rows.append({

                    "SMILES": smi, "Database": "SEA", "Rank": 0,

                    "Target_Name": "Processing failed", "Gene_Symbol": "",

                    "UniProt_ID": "", "ChEMBL_ID": "", "Target_Class": "",

                    "Target_Key": "", "Probability": "", "P_Value": "", "Max_Tc": "",

                })

            if idx < len(smiles_list):

                wait_time = random.uniform(10, 20)

                print(f"\nWaiting {wait_time:.2f} seconds before next compound...")

                time.sleep(wait_time)

    finally:

        driver.quit()

    # Save results with proper formatting

    columns = [

        "SMILES", "Database", "Target_Name", "Gene_Symbol",

        "UniProt_ID", "ChEMBL_ID", "Target_Class", "Target_Key",

        "Probability", "P_Value", "Max_Tc", "Rank"

    ]

    if swisstarget_rows:

        df_swiss = pd.DataFrame(swisstarget_rows)

        df_swiss = df_swiss.reindex(columns=columns)

        df_swiss = df_swiss.sort_values(["SMILES", "Probability"], ascending=[True, False])

        df_swiss['Rank'] = df_swiss.groupby('SMILES').cumcount() + 1

        df_swiss = df_swiss.sort_values(["SMILES", "Rank"]).reset_index(drop=True)

        df_swiss = df_swiss.fillna('')

        swiss_file = output_dir / "swisstarget_results.csv"

        df_swiss.to_csv(swiss_file, index=False)

        print(f"\n[SwissTarget] Saved {len(df_swiss)} rows to {swiss_file}")

    if sea_rows:

        df_sea = pd.DataFrame(sea_rows)

        df_sea = df_sea.reindex(columns=columns)

        df_sea = df_sea.sort_values(["SMILES", "P_Value"], ascending=[True, True])

        df_sea['Rank'] = df_sea.groupby('SMILES').cumcount() + 1

        df_sea = df_sea.sort_values(["SMILES", "Rank"]).reset_index(drop=True)

        df_sea = df_sea.fillna('')

        sea_file = output_dir / "sea_results.csv"

        df_sea.to_csv(sea_file, index=False)

        print(f"[SEA] Saved {len(df_sea)} rows to {sea_file}")

    # Combined results with proper formatting

    if swisstarget_rows or sea_rows:

        combined_file = output_dir / "combined_target_predictions.csv"

        format_combined_results(swisstarget_rows, sea_rows, combined_file)


def main():

    parser = argparse.ArgumentParser(

        description="Combined Target Prediction Tool (SwissTargetPrediction + SEA) - FORMATTING FIXED",

        formatter_class=argparse.RawDescriptionHelpFormatter

    )

    parser.add_argument(

        "input_file",

        help="Path to a text file containing one SMILES string per line"

    )

    parser.add_argument(

        "-o", "--output",

        default="results",

        help="Output directory for results (default: results)"

    )

    parser.add_argument(

        "--headless",

        action="store_true",

        help="Run browser in headless mode (no GUI)"

    )

    args = parser.parse_args()

    try:

        with open(args.input_file, 'r') as f:

            smiles_list = [line.strip() for line in f if line.strip()]

        if not smiles_list:

            print("Error: Input file is empty or contains no valid SMILES strings")

            return

    except FileNotFoundError:

        print(f"Error: Input file not found at {args.input_file}")

        return

    print(f"Loaded {len(smiles_list)} SMILES from {args.input_file}")

    print(f"Predictions will be filtered for HUMAN targets only (by UniProt ID suffix)\n")

    run_combined_prediction(smiles_list, output_dir=args.output, headless=args.headless)

    print("\n" + "="*60)

    print("Target prediction complete!")

    print("="*60)


if __name__ == "__main__":

    main()
