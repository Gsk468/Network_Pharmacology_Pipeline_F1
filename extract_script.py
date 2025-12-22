import json
import re

def extract_python_from_notebook(nb_path, out_path):
    with open(nb_path, "r", encoding="utf-8") as f:
        nb = json.load(f)

    script_content = []

    script_content.append("import os")
    script_content.append("import subprocess")
    script_content.append("import sys")
    script_content.append("# Mock get_ipython for !bang commands if we execute blindly, but better to replace them")
    script_content.append("\n")

    for i, cell in enumerate(nb["cells"]):
        if cell["cell_type"] == "code":
            source = cell["source"]
            source_str = "".join(source)

            # Heuristic to detect R cells
            if "library(" in source_str or "install.packages" in source_str or "ggsave" in source_str:
                script_content.append(f"# --- Cell {i}: Skipped (Appears to be R code) ---")
                continue

            script_content.append(f"# --- Cell {i} ---")

            for line in source:
                line_stripped = line.strip()
                if line_stripped.startswith("!"):
                    # Handle !pip install or !python
                    cmd = line_stripped[1:]
                    if "pip install" in cmd:
                        # Skip pip install in this runner, or warn
                        script_content.append(f"# Skipped: {line_stripped}")
                    elif "python" in cmd:
                        if "{" in cmd:
                             script_content.append(f"os.system(f'{cmd}')")
                        else:
                             script_content.append(f"os.system('{cmd}')")
                else:
                    script_content.append(line.rstrip())
            script_content.append("\n")

    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(script_content))
    print(f"Extracted Python script to {out_path}")

extract_python_from_notebook("Network_Pharmacology_Pipeline.ipynb", "pipeline_runner.py")
