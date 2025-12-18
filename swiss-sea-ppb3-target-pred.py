#!/usr/bin/env python3
"""
Target Prediction (Human-only) on 3 servers:
  1) SwissTargetPrediction (Homo sapiens)  [Selenium]
  2) SEA (Similarity Ensemble Approach)   [Selenium] -> filters by UniProt suffix "_HUMAN"
  3) PPB3                                 [requests] -> filters by Organism == Homo sapiens

Input:
  Sample-Plant-smiles.csv with columns:
    Phytochemical, Plant, SMILES, CID

Output (human targets only; failures are NOT written as rows):
  results_all3_human/
    combined.csv
    swiss.csv
    sea.csv
    ppb3.csv
    combined.json
"""

import os
import re
import csv
import time
import json
import random
import argparse
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Any, Optional

import pandas as pd

# --- Selenium (SwissTarget + SEA) ---
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from webdriver_manager.chrome import ChromeDriverManager

# --- PPB3 (requests) ---
import requests

try:
    from bs4 import BeautifulSoup
    HAS_BS4 = True
except ImportError:
    HAS_BS4 = False

# =========================
# UniProt Helpers
# =========================
UNIPROT_CACHE = {}

def fetch_uniprot_id(name: str) -> str:
    """
    Fetches UniProt ID (Entry Name or Accession) for a given gene/protein name.
    Target: Homo sapiens (9606).
    """
    if not name:
        return ""
    if name in UNIPROT_CACHE:
        return UNIPROT_CACHE[name]

    # Try searching as gene name first, then general query
    queries = [
        f"gene_exact:{name} AND organism_id:9606 AND reviewed:true",
        f"{name} AND organism_id:9606 AND reviewed:true",
        f"{name} AND organism_id:9606"
    ]

    api_url = "https://rest.uniprot.org/uniprotkb/search"

    for q in queries:
        try:
            params = {"query": q, "fields": "id,uniProtkbId", "size": 1, "format": "json"}
            resp = requests.get(api_url, params=params, timeout=10)
            if resp.status_code == 200:
                data = resp.json()
                results = data.get("results", [])
                if results:
                    # Prefer uniProtkbId (Entry Name like ALBU_HUMAN)
                    uid = results[0].get("uniProtkbId", "")
                    if not uid:
                        uid = results[0].get("primaryAccession", "")
                    if uid:
                        UNIPROT_CACHE[name] = uid
                        return uid
        except Exception:
            pass
        time.sleep(0.2)

    UNIPROT_CACHE[name] = "" # Cache failure to avoid retrying
    return ""



# =========================
# Configuration
# =========================
SWISSTARGET_URL = "https://www.swisstargetprediction.ch"
SEA_URL = "https://sea.bkslab.org/"
PPB3_URL = "https://ppb3.gdb.tools/"

USER_AGENTS = [
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/109.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/109.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/108.0.0.0 Safari/537.36",
    "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) AppleWebKit/605.1.15 (KHTML, like Gecko) Version/16.1 Safari/605.1.15",
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/108.0.0.0 Safari/537.36",
]

PPB3_PREDICTION_METHODS = [
    "DNN(ECFP4+MHFP6)",
    "DNN(ECFP4)",
    "DNN(RDKit)",
    "DNN(Layered)",
    "DNN(MHFP6)",
    "DNN(ECFP6)",
    "DNN(AtomPair)",
    "Consensus",
]

PPB3_MAX_RETRIES = 3
PPB3_TIMEOUT = 60


# =========================
# Helpers: human filters
# =========================
def is_human_target_uniprot(target_key: str) -> bool:
    """
    SEA human filtering:
    Human targets end with "_HUMAN" in Target_Key (UniProt-style ID).
    """
    if not target_key or not isinstance(target_key, str):
        return False
    return target_key.upper().endswith("_HUMAN")


def is_homo_sapiens(text: str) -> bool:
    """
    PPB3 organism filter:
    Keep only Homo sapiens / H. sapiens.
    """
    if not text or not isinstance(text, str):
        return False
    t = text.strip().lower()
    return (t == "homo sapiens") or (t == "h. sapiens") or ("homo sapiens" in t)


# =========================
# Input
# =========================
def read_phytochemical_csv(csv_path: str, max_compounds: Optional[int] = None) -> List[Dict[str, str]]:
    """
    Reads Sample-Plant-smiles.csv-like file.
    Requires: Phytochemical, SMILES
    Optional: Plant, CID
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    out = []
    with open(csv_path, "r", encoding="utf-8", errors="ignore") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames or []
        lower_map = {c.lower().strip(): c for c in fieldnames}

        phy_col = lower_map.get("phytochemical")
        plant_col = lower_map.get("plant")
        smi_col = lower_map.get("smiles")
        cid_col = lower_map.get("cid")

        if not phy_col or not smi_col:
            raise RuntimeError(f"Missing required columns (need Phytochemical + SMILES). Found: {fieldnames}")

        for row in reader:
            if max_compounds and len(out) >= max_compounds:
                break
            phy = (row.get(phy_col) or "").strip()
            smi = (row.get(smi_col) or "").strip()
            if not phy or not smi:
                continue

            out.append({
                "Phytochemical": phy,
                "Plant": (row.get(plant_col) or "").strip() if plant_col else "",
                "SMILES": smi,
                "CID": (row.get(cid_col) or "").strip() if cid_col else "",
            })

    return out


# =========================
# Selenium setup
# =========================
def setup_driver(headless: bool = False) -> webdriver.Chrome:
    opts = Options()
    ua = random.choice(USER_AGENTS)
    opts.add_argument(f"user-agent={ua}")
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
    driver.execute_script("Object.defineProperty(navigator, 'webdriver', {get: () => undefined})")
    return driver


def human_like_typing(element, text: str):
    for ch in text:
        element.send_keys(ch)
        time.sleep(random.uniform(0.05, 0.2))


# =========================
# SwissTargetPrediction (Homo sapiens)
# =========================
def swisstarget_submit(driver, smiles: str, species: str = "Homo sapiens"):
    driver.get(SWISSTARGET_URL + "/")
    WebDriverWait(driver, 30).until(EC.presence_of_element_located((By.TAG_NAME, "body")))
    time.sleep(random.uniform(1, 3))

    smiles_input = WebDriverWait(driver, 20).until(
        EC.presence_of_element_located((By.CSS_SELECTOR, "input[name='smiles']"))
    )
    smiles_input.clear()
    human_like_typing(smiles_input, smiles)
    time.sleep(random.uniform(0.5, 1.5))

    # Select organism/species (human)
    try:
        species_sel = driver.find_element(By.CSS_SELECTOR, "select[name='species'], select[name='organism']")
        species_sel.send_keys(species)
        time.sleep(random.uniform(0.5, 1.5))
    except Exception:
        pass

    submit_button_selectors = ["input[type='submit']", "input[value*='Predict']", "#submitbutton"]
    submit = None
    for sel in submit_button_selectors:
        try:
            submit = driver.find_element(By.CSS_SELECTOR, sel)
            break
        except Exception:
            continue
    if not submit:
        raise RuntimeError("[SwissTarget] Could not find submit button")

    driver.execute_script("arguments[0].click();", submit)


def swisstarget_get_results(driver, max_wait: int = 300):
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


def swisstarget_parse_table(table, meta: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    SwissTarget is already set to Homo sapiens by the submission species selector.
    Parses the results table into standard rows.
    """
    smiles = meta["SMILES"]
    rows = []

    headers = [th.text.strip().lower() for th in table.find_elements(By.TAG_NAME, "th")]
    header_map = {
        "target": next((h for h in headers if "target" in h), None),
        "gene": next((h for h in headers if "common" in h), None),
        "uniprot": next((h for h in headers if "uniprot" in h), None),
        "chembl": next((h for h in headers if "chembl" in h), None),
        "class": next((h for h in headers if "class" in h), None),
        "prob": next((h for h in headers if "probability" in h), None),
    }
    col_indices = {k: headers.index(v) for k, v in header_map.items() if v is not None}

    prob_regex = re.compile(r"[0-9]*\.?[0-9]+")
    trs = table.find_elements(By.TAG_NAME, "tr")

    for i, tr in enumerate(trs[1:], start=1):
        tds = tr.find_elements(By.TAG_NAME, "td")
        if len(tds) < len(col_indices):
            continue

        def get_text(col):
            idx = col_indices.get(col)
            return tds[idx].text.strip() if idx is not None and idx < len(tds) else ""

        target_name = get_text("target")
        uniprot_id = get_text("uniprot")
        if not target_name or not uniprot_id:
            continue

        prob_text = get_text("prob")
        m = prob_regex.search(prob_text)
        prob_val = float(m.group(0)) if m else 0.0

        rows.append({
            "Phytochemical": meta["Phytochemical"],
            "Plant": meta["Plant"],
            "CID": meta["CID"],
            "SMILES": smiles,
            "Database": "SwissTargetPrediction",
            "Rank": i,
            "Target_Name": target_name,
            "Gene_Symbol": get_text("gene"),
            "UniProt_ID": uniprot_id,
            "ChEMBL_ID": get_text("chembl"),
            "Target_Class": get_text("class"),
            "Target_Key": uniprot_id,
            "Probability": prob_val,
            "P_Value": "",
            "Max_Tc": "",
        })

    return rows


# =========================
# SEA (Human only by _HUMAN suffix)
# =========================
def sea_submit(driver, smiles: str):
    driver.get(SEA_URL + "/")
    WebDriverWait(driver, 30).until(EC.presence_of_element_located((By.TAG_NAME, "body")))
    time.sleep(random.uniform(1, 3))

    smiles_input = WebDriverWait(driver, 20).until(
        EC.presence_of_element_located((By.CSS_SELECTOR, "input[placeholder='Paste SMILES or try the example below']"))
    )
    smiles_input.clear()
    human_like_typing(smiles_input, smiles)
    time.sleep(random.uniform(0.5, 1.5))

    submit_button = WebDriverWait(driver, 10).until(
        EC.element_to_be_clickable((By.XPATH, "//button[contains(., 'Try SEA')]"))
    )
    submit_button.click()


def sea_get_results(driver, max_wait: int = 300):
    WebDriverWait(driver, max_wait).until(EC.url_contains("/jobs/"))
    table = WebDriverWait(driver, 60).until(
        EC.presence_of_element_located((By.CSS_SELECTOR, "table.table.table-bordered"))
    )
    return table


def sea_parse_table(table, meta: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Parses SEA results and keeps only human targets (Target_Key ends with _HUMAN).
    """
    smiles = meta["SMILES"]
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
    for out_col, header_text in col_mapping.items():
        idx = thead_header_texts.index(header_text) - 1  # SEA has a leading blank/rank column
        actual_td_indices[out_col] = idx

    tbody_rows = table.find_elements(By.TAG_NAME, "tbody")[0].find_elements(By.TAG_NAME, "tr")

    rank = 0
    for tr in tbody_rows:
        tr_class = tr.get_attribute("class")
        if tr_class and "spanning info" in tr_class:
            continue

        tds = tr.find_elements(By.TAG_NAME, "td")
        if len(tds) < len(col_mapping):
            continue

        target_key = tds[actual_td_indices["Target_Key"]].text.strip()
        if not is_human_target_uniprot(target_key):
            continue

        rank += 1
        target_name = tds[actual_td_indices["Target_Name"]].text.strip()
        desc = tds[actual_td_indices["Description"]].text.strip()

        rows.append({
            "Phytochemical": meta["Phytochemical"],
            "Plant": meta["Plant"],
            "CID": meta["CID"],
            "SMILES": smiles,
            "Database": "SEA",
            "Rank": rank,
            "Target_Name": target_name,
            "Gene_Symbol": "",
            "UniProt_ID": target_key,
            "ChEMBL_ID": "",
            "Target_Class": desc,
            "Target_Key": target_key,
            "Probability": float(tds[actual_td_indices["P_Value"]].text.strip()) if tds[actual_td_indices["P_Value"]].text.strip() else 1.0, # Store P-value in Probability column
            "P_Value": tds[actual_td_indices["P_Value"]].text.strip(),
            "Max_Tc": tds[actual_td_indices["Max_Tc"]].text.strip(),
        })

    return rows


# =========================
# PPB3 (requests + HTML parse) + Homo sapiens only
# =========================
def ppb3_extract_targets_from_page(page_html: str) -> List[Dict[str, Any]]:
    """
    Returns list of targets parsed from PPB3 HTML response.
    Human-only filter applied based on Organism column (cell[5]).
    """
    out: List[Dict[str, Any]] = []
    if not HAS_BS4:
        return out

    soup = BeautifulSoup(page_html, "html.parser")
    table = soup.find("table", id="resultsTable")
    if not table:
        return out

    for row in table.select("tbody tr"):
        cells = row.find_all("td")
        if len(cells) < 6:
            continue

        try:
            target_name = cells[2].get_text(strip=True)
            prob_str = cells[3].get_text(strip=True)
            organism = cells[5].get_text(strip=True)

            # HUMAN ONLY (Homo sapiens)
            if not is_homo_sapiens(organism):
                continue

            prob = float(prob_str)
            if target_name:
                out.append({"Target_Name": target_name, "Probability": prob, "Organism": organism})
        except Exception:
            continue

    return out


def ppb3_submit_smiles_for_method(smiles: str, method: str, retry: int = 0) -> List[Dict[str, Any]]:
    if retry >= PPB3_MAX_RETRIES:
        return []

    model_type_map = {
        "DNN(ECFP4+MHFP6)": "Fused",
        "DNN(ECFP4)": "ECFP4",
        "DNN(RDKit)": "RDKit",
        "DNN(Layered)": "Layered",
        "DNN(MHFP6)": "MHFP6",
        "DNN(ECFP6)": "ECFP6",
        "DNN(AtomPair)": "AtomPair",
        "Consensus": "Consensus",
    }
    model_type = model_type_map.get(method)
    if not model_type:
        return []

    url = f"{PPB3_URL}result"
    payload = {"smiles": [smiles], "model_type": model_type}

    try:
        r = requests.post(url, json=payload, timeout=PPB3_TIMEOUT)
        r.raise_for_status()
        return ppb3_extract_targets_from_page(r.text)
    except Exception:
        time.sleep(2 * (retry + 1))
        return ppb3_submit_smiles_for_method(smiles, method, retry + 1)


def run_ppb3_for_compound(meta: Dict[str, str]) -> List[Dict[str, Any]]:
    """
    Runs PPB3 across all methods and returns only Homo sapiens targets.
    Notes:
      - Target_Class stores PPB3 method name.
      - Target_Key stores organism string.
    """
    smiles = meta["SMILES"]
    rows: List[Dict[str, Any]] = []

    for method in PPB3_PREDICTION_METHODS:
        targets = ppb3_submit_smiles_for_method(smiles, method)
        for t in targets:
            rows.append({
                "Phytochemical": meta["Phytochemical"],
                "Plant": meta["Plant"],
                "CID": meta["CID"],
                "SMILES": smiles,
                "Database": "PPB3",
                "Rank": "",  # recalculated later
                "Target_Name": t.get("Target_Name", ""),
                "Gene_Symbol": "",
                "UniProt_ID": fetch_uniprot_id(t.get("Target_Name", "")),
                "ChEMBL_ID": "",
                "Target_Class": method,                 # PPB3 model/method
                "Target_Key": t.get("Organism", ""),    # organism
                "Probability": t.get("Probability", ""),
                "P_Value": "",
                "Max_Tc": "",
            })

    return rows


# =========================
# Output formatting
# =========================
OUTPUT_COLUMNS = [
    "Phytochemical", "Plant", "CID", "SMILES",
    "Database", "Target_Name", "Gene_Symbol", "UniProt_ID", "ChEMBL_ID",
    "Target_Class", "Target_Key", "Probability", "P_Value", "Max_Tc", "Rank"
]


def finalize_rank(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["Rank"] = df.groupby(["Phytochemical", "SMILES", "Database"]).cumcount() + 1
    return df


def save_outputs(all_rows: List[Dict[str, Any]], output_dir: str):
    outdir = Path(output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = pd.DataFrame(all_rows)
    if df.empty:
        print("No HUMAN targets collected from any server; nothing to write.")
        return

    df = df.reindex(columns=OUTPUT_COLUMNS).fillna("")
    df = finalize_rank(df)
    df = df.sort_values(["Phytochemical", "Database", "Rank"]).reset_index(drop=True)

    combined_csv = outdir / "combined.csv"
    df.to_csv(combined_csv, index=False)
    print(f"[OK] Combined HUMAN CSV: {combined_csv} ({len(df)} rows)")

    for db in sorted(df["Database"].unique()):
        sub = df[df["Database"] == db].copy()
        out_csv = outdir / f"{db.lower()}_results_human.csv"
        sub.to_csv(out_csv, index=False)
        print(f"[OK] {db} HUMAN CSV: {out_csv} ({len(sub)} rows)")

    combined_json = outdir / "combined.json"
    with open(combined_json, "w", encoding="utf-8") as f:
        json.dump(
            {
                "metadata": {
                    "timestamp": datetime.now().isoformat(),
                    "input": "CSV (Phytochemical, Plant, SMILES, CID)",
                    "databases": ["SwissTargetPrediction", "SEA", "PPB3"],
                    "human_filtering": {
                        "SwissTargetPrediction": "Submitted with species=Homo sapiens",
                        "SEA": "Target_Key endswith _HUMAN",
                        "PPB3": "Organism contains Homo sapiens",
                    },
                },
                "rows": df.to_dict(orient="records"),
            },
            f,
            indent=2,
        )
    print(f"[OK] Combined HUMAN JSON: {combined_json}")


# =========================
# Main workflow
# =========================
def main():
    ap = argparse.ArgumentParser(description="Human-only target prediction on SwissTargetPrediction + SEA + PPB3 from CSV.")
    ap.add_argument("csv_file", help="Path to Sample-Plant-smiles.csv")
    ap.add_argument("-o", "--output", default="results_all3_human", help="Output directory (default: results_all3_human)")
    ap.add_argument("--headless", action="store_true", help="Run Selenium in headless mode")
    ap.add_argument("--max-compounds", type=int, default=None, help="Process only first N compounds")
    ap.add_argument("--min-wait", type=float, default=10.0, help="Min seconds between compounds (Swiss/SEA)")
    ap.add_argument("--max-wait", type=float, default=20.0, help="Max seconds between compounds (Swiss/SEA)")
    args = ap.parse_args()

    compounds = read_phytochemical_csv(args.csv_file, max_compounds=args.max_compounds)
    if not compounds:
        print("No compounds found in CSV (or missing SMILES).")
        return

    if not HAS_BS4:
        print("WARNING: beautifulsoup4 not installed; PPB3 parsing will return 0 rows. Install with: pip install beautifulsoup4")

    print(f"Loaded {len(compounds)} compounds from {args.csv_file}")
    print("Filtering to HUMAN (Homo sapiens) targets only for all 3 servers.\n")

    all_rows: List[Dict[str, Any]] = []

    # Try setting up driver, proceed without it if it fails (run only PPB3)
    try:
        driver = setup_driver(headless=args.headless)
    except Exception as e:
        print(f"WARNING: Could not set up Selenium Driver (Chrome). Skipping Swiss/SEA. Error: {e}")
        driver = None

    try:
        for i, meta in enumerate(compounds, start=1):
            print("\n" + "=" * 70)
            print(f"[{i}/{len(compounds)}] {meta['Phytochemical']} | {meta['Plant']}")
            print(f"SMILES: {meta['SMILES']}", flush=True)
            print("=" * 70)

            # SwissTargetPrediction (already Homo sapiens)
            if driver:
                try:
                    print("[SwissTarget] Running (Homo sapiens)...", flush=True)
                    swisstarget_submit(driver, meta["SMILES"], species="Homo sapiens")
                    table = swisstarget_get_results(driver)
                    rows = swisstarget_parse_table(table, meta)
                    print(f"[SwissTarget] Parsed {len(rows)} human rows", flush=True)
                    all_rows.extend(rows)
                except Exception as e:
                    print(f"[SwissTarget] FAILED (skipping rows): {e}", flush=True)

                time.sleep(random.uniform(3, 6))

                # SEA (human only by _HUMAN)
                try:
                    print("[SEA] Running (human-only by _HUMAN)...", flush=True)
                    sea_submit(driver, meta["SMILES"])
                    table = sea_get_results(driver)
                    rows = sea_parse_table(table, meta)
                    print(f"[SEA] Parsed {len(rows)} human rows", flush=True)
                    all_rows.extend(rows)
                except Exception as e:
                    print(f"[SEA] FAILED (skipping rows): {e}", flush=True)
            else:
                print("[SwissTarget/SEA] Skipped (No Selenium Driver)", flush=True)

            # PPB3 (human only by Organism)
            try:
                print("[PPB3] Running (Homo sapiens only by Organism column)...", flush=True)
                rows = run_ppb3_for_compound(meta)
                print(f"[PPB3] Parsed {len(rows)} human rows (across methods)", flush=True)
                all_rows.extend(rows)
            except Exception as e:
                print(f"[PPB3] FAILED (skipping rows): {e}", flush=True)

            if i < len(compounds):
                wait = random.uniform(args.min_wait, args.max_wait)
                print(f"Waiting {wait:.2f}s before next compound...", flush=True)
                time.sleep(wait)

    finally:
        if driver:
            driver.quit()

    save_outputs(all_rows, args.output)
    print("Done.")


if __name__ == "__main__":
    main()
