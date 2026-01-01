import argparse
import pandas as pd
import os
import random

# List of real human gene symbols (Uniprot ID | Symbol format for targets.py)
REAL_TARGETS = [
    'TP53_HUMAN|P04637', 'TNF_HUMAN|P01375', 'EGFR_HUMAN|P00533',
    'AKT1_HUMAN|P31749', 'IL6_HUMAN|P05231', 'VEGFA_HUMAN|P15692',
    'GAPDH_HUMAN|P04406', 'INS_HUMAN|P01308', 'ALB_HUMAN|P02768',
    'CTNNB1_HUMAN|P35222', 'MYC_HUMAN|P01106', 'JUN_HUMAN|P05412'
]

def main():
    parser = argparse.ArgumentParser(description='Mock targets.py for testing')
    parser.add_argument('--input', required=True, help='Input CSV file')
    parser.add_argument('--output', required=True, help='Output CSV file')
    args = parser.parse_args()

    print(f'Mock targets.py: Processing {args.input} to {args.output}')

    if os.path.exists(args.input):
        try:
            df_in = pd.read_csv(args.input)
            results = []
            for _, row in df_in.iterrows():
                ing = row.get('Ingredient', 'Unknown')
                selected_targets = random.sample(REAL_TARGETS, k=min(5, len(REAL_TARGETS)))
                for target in selected_targets:
                    results.append({
                        'Ingredient': ing,
                        'UniProt_name': target,
                        'Probability': round(random.uniform(0.5, 1.0), 2)
                    })
            df_out = pd.DataFrame(results)
            os.makedirs(os.path.dirname(args.output), exist_ok=True)
            df_out.to_csv(args.output, index=False)
            print(f'Mock targets.py: Generated {len(df_out)} targets.')
        except Exception as e:
            print(f'Error in mock targets.py: {e}')
    else:
        print(f'Error: Input file {args.input} not found.')

if __name__ == '__main__':
    main()