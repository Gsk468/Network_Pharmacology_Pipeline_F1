#!/usr/bin/env python3
import os,sys,re,csv,json,time,logging
from pathlib import Path
from typing import Dict,Set,Optional
from collections import defaultdict
from datetime import datetime
import argparse,traceback

# --- New dependency ---
try:
    import requests
except ImportError:
    print("The 'requests' library is not installed. Please install it by running: pip install requests")
    sys.exit(1)

logging.basicConfig(level=logging.INFO,format='%(asctime)s - %(levelname)s - %(message)s',handlers=[logging.FileHandler('ppb3_automation.log'),logging.StreamHandler()])
logger = logging.getLogger(__name__)

try:
    from bs4 import BeautifulSoup;HAS_BS4 = True
except ImportError:
    HAS_BS4 = False

PPB3_URL = "https://ppb3.gdb.tools/"
PREDICTION_METHODS = ["DNN(ECFP4+MHFP6)","DNN(ECFP4)","DNN(RDKit)","DNN(Layered)","DNN(MHFP6)","DNN(ECFP6)","DNN(AtomPair)","Consensus"]
MAX_RETRIES,TIMEOUT,INTER_COMPOUND_DELAY,INTER_METHOD_DELAY = 3,60,1,1 # Adjusted delays for API calls

class PPB3Automation:
    def __init__(self,csv_path: str,max_compounds: Optional[int]=None,output_dir: str="ppb3_results",headless: bool=False): # headless is no longer used but kept for CLI compatibility
        self.csv_path = csv_path
        self.max_compounds = max_compounds
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True,parents=True)
        self.smiles_data = {}
        self.results = defaultdict(lambda: defaultdict(set))
        self.failed_smiles = []
        self.start_time = None

    def parse_csv(self):
        smiles_dict, seen_smiles = {}, set()
        logger.info(f"Parsing CSV: {self.csv_path}")
        if not os.path.exists(self.csv_path):
            logger.error(f"CSV file not found: {self.csv_path}")
            return False
        try:
            with open(self.csv_path, 'r', encoding='utf-8', errors='ignore') as f:
                reader = csv.DictReader(f)
                headers = [h.lower().strip() for h in reader.fieldnames]
                name_headers = ['name', 'compound', 'compoundname', 'id', 'compoundid']
                smiles_headers = ['smiles', 'smi', 'canonicalsmiles', 'canonical smiles']
                
                name_col = next((reader.fieldnames[headers.index(h)] for h in name_headers if h in headers), None)
                smiles_col = next((reader.fieldnames[headers.index(h)] for h in smiles_headers if h in headers), None)

                if not name_col or not smiles_col:
                    logger.error(f"Could not find a 'Name' and/or 'SMILES' column in the CSV header. Found: {reader.fieldnames}")
                    return False
                
                logger.info(f"Using column '{name_col}' for names and '{smiles_col}' for SMILES.")
                
                rows = list(reader)
                total_lines = len(rows)
                logger.info(f"Found {total_lines} data rows to process.")
                
                count = 0
                for row in rows:
                    if self.max_compounds and count >= self.max_compounds: 
                        logger.info(f"Reached max compounds limit: {self.max_compounds}"); break
                    
                    name = row[name_col].strip()
                    smiles = row[smiles_col].strip()

                    if not name or not smiles: continue
                    if not self.is_valid_smiles(smiles): continue
                    if smiles in seen_smiles: continue
                        
                    smiles_dict[name] = smiles
                    seen_smiles.add(smiles)
                    count += 1
            
            logger.info(f"Successfully parsed {len(smiles_dict)} unique compounds.")
            self.smiles_data = smiles_dict
            return True
        except Exception as e:
            logger.error(f"Fatal error parsing CSV: {e}")
            logger.error(traceback.format_exc())
            return False

    def is_valid_smiles(self,smiles: str)->bool:
        if not smiles or len(smiles)<2:return False
        valid_elements = set('CBNOSPFClBrI')
        return any(elem in smiles for elem in valid_elements)

    def extract_targets_from_page(self, page_html: str) -> Set[tuple[str, float, str]]:
        # Returns a set of tuples: (target_name, probability, organism)
        targets = set()
        try:
            if HAS_BS4:
                soup = BeautifulSoup(page_html, 'html.parser')
                table = soup.find('table', id='resultsTable')
                if table:
                    for row in table.select('tbody tr'):
                        cells = row.find_all('td')
                        if len(cells) >= 6: # Need at least 6 cells for all data
                            try:
                                target_name = cells[2].get_text(strip=True)
                                probability_str = cells[3].get_text(strip=True)
                                organism = cells[5].get_text(strip=True)
                                
                                if target_name and probability_str:
                                    probability = float(probability_str)
                                    targets.add((target_name, probability, organism))
                            except (ValueError, IndexError):
                                logger.warning("Could not parse a row in the results table.")
                                continue
                else:
                    logger.warning("Could not find results table in response HTML. Prediction may have failed on the server.")
            else:
                logger.warning("BeautifulSoup not installed. Cannot parse HTML results.")
        except Exception as e:
            logger.error(f"Error parsing results HTML: {e}")
        return targets

    def submit_smiles_for_method(self,smiles: str,method: str,retry: int=0)->Set[str]:
        if retry >= MAX_RETRIES:
            logger.error(f"Max retries exceeded for {smiles[:30]}... with method {method}"); return set()

        model_type_map = {
            "DNN(ECFP4+MHFP6)": "Fused", "DNN(ECFP4)": "ECFP4", "DNN(RDKit)": "RDKit",
            "DNN(Layered)": "Layered", "DNN(MHFP6)": "MHFP6", "DNN(ECFP6)": "ECFP6",
            "DNN(AtomPair)": "AtomPair", "Consensus": "Consensus"
        }
        model_type = model_type_map.get(method)
        if not model_type:
            logger.error(f"Unknown method: {method}"); return set()

        url = f"{PPB3_URL}result"
        payload = {"smiles": [smiles], "model_type": model_type}
        
        try:
            response = requests.post(url, json=payload, timeout=TIMEOUT)
            response.raise_for_status()
            return self.extract_targets_from_page(response.text)
        except requests.exceptions.RequestException as e:
            logger.warning(f"HTTP Request failed for method {method} (attempt {retry+1}/{MAX_RETRIES}): {e}")
            time.sleep(2 * (retry + 1)) # Exponential backoff
            return self.submit_smiles_for_method(smiles, method, retry + 1)

    def run_predictions(self):
        total = len(self.smiles_data)
        if total == 0:
            logger.error("No compounds to process!");return False
        self.start_time = datetime.now()
        
        try:
            for idx,(compound_name,smiles) in enumerate(self.smiles_data.items(),1):
                logger.info("="*70); logger.info(f"Processing compound {idx}/{total}: {compound_name}"); logger.info(f"SMILES: {smiles}")
                
                compound_success = False
                for method in PREDICTION_METHODS:
                    logger.info(f"Running method '{method}'...")
                    targets = self.submit_smiles_for_method(smiles,method)
                    if targets:
                        self.results[compound_name][method] = targets
                        compound_success = True
                    else:
                        logger.warning(f"Method '{method}' returned no targets.")
                    time.sleep(INTER_METHOD_DELAY)
                
                if not compound_success:
                    self.failed_smiles.append((compound_name,smiles))
                    logger.warning(f"Failed to get any results for: {compound_name}")
                
                logger.info(f"Finished processing compound {idx}/{total}: {compound_name}")
                time.sleep(INTER_COMPOUND_DELAY)
        except KeyboardInterrupt:
            logger.info("Interrupted by user")
        return True

    def generate_consensus(self):
        consensus = {}
        for compound_name, method_results in self.results.items():
            all_target_names = set()
            for targets_with_details in method_results.values():
                for target_info in targets_with_details:
                    all_target_names.add(target_info[0]) # Add just the name
            consensus[compound_name] = all_target_names
        return consensus

    def save_results(self):
        logger.info("="*70); logger.info("SAVING RESULTS"); logger.info("="*70)
        consensus = self.generate_consensus()
        json_file = self.output_dir / "ppb3_predictions.json"
        try:
            all_targets = set()
            for targets in consensus.values(): all_targets.update(targets)
            avg_targets = sum(len(t) for t in consensus.values()) / len(consensus) if consensus else 0
            report = {
                "metadata": {
                    "timestamp": datetime.now().isoformat(), "source": "PPB3", "url": PPB3_URL,
                    "methods": PREDICTION_METHODS, "total_compounds": len(self.smiles_data),
                    "successful": len(self.results), "failed": len(self.failed_smiles)
                },
                "summary": {
                    "unique_targets": len(all_targets), "avg_targets_per_compound": round(avg_targets, 2),
                    "max_targets": max((len(t) for t in consensus.values()), default=0),
                    "min_targets": min((len(t) for t in consensus.values()), default=0)
                },
                "results": {}
            }
            for compound in self.smiles_data:
                report["results"][compound] = {
                    "smiles": self.smiles_data[compound],
                    "targets_by_method": {method: sorted([f"{t[0]} ({t[1]:.2f}, {t[2]})" for t in targets]) for method, targets in self.results.get(compound, {}).items()},
                    "consensus_targets": sorted(list(consensus.get(compound, set())))
                }
            with open(json_file, 'w') as f: json.dump(report, f, indent=2)
            logger.info(f"JSON report saved: {json_file}")
        except Exception as e:
            logger.error(f"Error saving JSON report: {e}")

        csv_file = self.output_dir / "ppb3_consensus_targets.csv"
        try:
            with open(csv_file, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["Compound Name", "SMILES", "Number of Targets", "Consensus Targets"])
                for compound_name, smiles in self.smiles_data.items():
                    targets = consensus.get(compound_name, set())
                    writer.writerow([compound_name, smiles, len(targets), "; ".join(sorted(targets))])
            logger.info(f"Consensus CSV saved: {csv_file}")
        except Exception as e:
            logger.error(f"Error saving consensus CSV: {e}")

        detailed_csv = self.output_dir / "ppb3_detailed_results.csv"
        try:
            with open(detailed_csv, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["Compound Name", "SMILES", "Method", "Targets (Probability, Organism)"])
                for compound_name, methods in self.results.items():
                    smiles = self.smiles_data.get(compound_name, "")
                    for method, targets in methods.items():
                        formatted_targets = [f"{t[0]} ({t[1]:.2f}, {t[2]})" for t in targets]
                        writer.writerow([compound_name, smiles, method, "; ".join(sorted(formatted_targets))])
            logger.info(f"Detailed CSV saved: {detailed_csv}")
        except Exception as e:
            logger.error(f"Error saving detailed CSV: {e}")

        # --- New User-Requested Output File ---
        full_details_csv = self.output_dir / "ppb3_full_details.csv"
        try:
            with open(full_details_csv, 'w', newline='', encoding='utf-8') as f:
                writer = csv.writer(f)
                writer.writerow(["Compound Name", "SMILES", "Method", "Target Name", "Organism", "Probability"])
                for compound_name, methods_results in self.results.items():
                    smiles = self.smiles_data.get(compound_name, "")
                    for method, targets in methods_results.items():
                        for target_info in targets:
                            target_name, probability, organism = target_info
                            writer.writerow([compound_name, smiles, method, target_name, organism, probability])
            logger.info(f"Full detailed CSV saved: {full_details_csv}")
        except Exception as e:
            logger.error(f"Error saving full detailed CSV: {e}")

        if self.failed_smiles:
            failed_csv = self.output_dir / "ppb3_failed_compounds.csv"
            try:
                with open(failed_csv, 'w', newline='', encoding='utf-8') as f:
                    writer = csv.writer(f)
                    writer.writerow(["Compound Name", "SMILES", "Reason"])
                    for name, smiles in self.failed_smiles:
                        writer.writerow([name, smiles, "Failed to retrieve predictions"])
                logger.info(f"Failed compounds CSV saved: {failed_csv}")
            except Exception as e:
                logger.error(f"Error saving failed compounds CSV: {e}")

        logger.info("="*70); logger.info(f"Total Compounds Processed: {len(self.smiles_data)}"); logger.info(f"Successful Predictions: {len(self.results)}"); logger.info(f"Failed Predictions: {len(self.failed_smiles)}"); logger.info("="*70)

    def print_summary(self):
        consensus = self.generate_consensus()
        print("\n" + "="*70); print("PPB3 BACTERIAL TARGET PREDICTION - SUMMARY"); print("="*70)
        print(f"Total Compounds: {len(self.smiles_data)}"); print(f"Successful: {len(self.results)}"); print(f"Failed: {len(self.failed_smiles)}")
        if self.start_time: print(f"Elapsed Time: {datetime.now() - self.start_time}")
        print(f"Methods Used: {len(PREDICTION_METHODS)}")
        for i, method in enumerate(PREDICTION_METHODS, 1): print(f"  {i}. {method}")
        if consensus:
            all_targets = set()
            for targets in consensus.values(): all_targets.update(targets)
            print(f"\nStatistics:"); print(f"  Total Unique Targets: {len(all_targets)}")
            target_counts = [len(t) for t in consensus.values()]
            avg = sum(target_counts) / len(target_counts) if target_counts else 0
            print(f"  Average Targets per Compound: {avg:.2f}"); print(f"  Max Targets: {max(target_counts) if target_counts else 0}"); print(f"  Min Targets: {min(target_counts) if target_counts else 0}")
        print(f"\nFiles Generated:"); print(f"  Location: {self.output_dir.absolute()}"); print(f"  - ppb3_predictions.json"); print(f"  - ppb3_consensus_targets.csv"); print(f"  - ppb3_detailed_results.csv")
        if self.failed_smiles: print(f"  - ppb3_failed_compounds.csv")
        print("="*70); print()

def main():
    parser = argparse.ArgumentParser(description='PPB3 Automated Bacterial Target Prediction', formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('csv_file', help='Path to CSV file containing SMILES data')
    parser.add_argument('--max-compounds', type=int, default=None, help='Maximum number of compounds to process')
    parser.add_argument('--output-dir', default='ppb3_results', help='Output directory for results')
    parser.add_argument('--headless', action='store_true', help='(No longer used) Kept for compatibility.')
    parser.add_argument('--debug', action='store_true', help='Enable debug logging')
    args = parser.parse_args()
    if args.debug:
        logger.setLevel(logging.DEBUG)
    
    automation = PPB3Automation(csv_path=args.csv_file, max_compounds=args.max_compounds, output_dir=args.output_dir)
    if not automation.parse_csv():
        sys.exit(1)
    
    automation.run_predictions()
    automation.save_results()
    automation.print_summary()
    logger.info("PPB3 AUTOMATION COMPLETED SUCCESSFULLY")
    sys.exit(0)

if __name__ == '__main__':
    main()