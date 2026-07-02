#!/usr/bin/env python3
"""
generate_citations_table.py
===========================
Extracts the citation data from the projects.csv file in every camera trap dataset
under the data/ directory, de-duplicates by project_id, and creates a neat,
publication-quality markdown table citing all data used.
"""

import os
import pandas as pd
from pathlib import Path

def make_markdown_table(df):
    headers = list(df.columns)
    lines = []
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for _, row in df.iterrows():
        val_str = []
        for h in headers:
            val = str(row[h]).replace("|", "\\|").replace("\n", " ").replace("\r", " ").strip()
            val_str.append(val)
        lines.append("| " + " | ".join(val_str) + " |")
    return "\n".join(lines)

def main():
    print("=== Extracting Camera Trap Citations ===")
    
    # 1. Find all projects.csv files recursively in data/
    data_dir = Path("data")
    projects_files = list(data_dir.glob("**/projects.csv"))
    print(f"Found {len(projects_files)} projects.csv files in data directory.")
    
    records = []
    
    # 2. Extract project info and citations
    for path in projects_files:
        try:
            # We use keep_default_na=False to avoid interpreting NA or None as NaN
            df = pd.read_csv(path, keep_default_na=False)
            
            # Check if required columns exist
            required = ["project_id", "project_name", "data_citation"]
            if not all(col in df.columns for col in required):
                # Fallback to look for case insensitive matches
                col_map = {col.lower(): col for col in df.columns}
                if all(req.lower() in col_map for req in required):
                    df = df.rename(columns={col_map[req.lower()]: req for req in required})
                else:
                    print(f"Warning: projects.csv at {path} is missing required columns. Skipping.")
                    continue
            
            for _, row in df.iterrows():
                project_id = str(row["project_id"]).strip()
                project_name = str(row["project_name"]).strip()
                data_citation = str(row["data_citation"]).strip()
                
                # Determine basin based on path
                path_str = str(path)
                if "Amazon" in path_str:
                    basin = "Amazon"
                elif "Congo" in path_str:
                    basin = "Congo"
                elif "SEAsia" in path_str:
                    basin = "SE Asia"
                else:
                    basin = "Other"
                
                if project_id and data_citation:
                    records.append({
                        "Project ID": project_id,
                        "Project Name": project_name,
                        "Basin": basin,
                        "Citation": data_citation
                    })
        except Exception as e:
            print(f"Error reading {path}: {e}")
            
    if not records:
        print("No citations found! Exiting.")
        return
        
    df_all = pd.DataFrame(records)
    
    # 3. De-duplicate by Project ID
    df_unique = df_all.drop_duplicates(subset=["Project ID"])
    df_unique = df_unique.sort_values(by=["Basin", "Project Name"]).reset_index(drop=True)
    
    print(f"Extracted {len(df_unique)} unique project citations out of {len(df_all)} total entries.")
    
    # 4. Generate Markdown Table
    df_md = df_unique[["Basin", "Project Name", "Project ID", "Citation"]].copy()
    
    md_content = "# Data Citations Table\n\n"
    md_content += "This table lists the unique camera trap projects used in the analysis, categorized by region, along with their respective metadata license and dataset citation.\n\n"
    md_content += make_markdown_table(df_md)
    md_content += "\n"
    
    # Write outputs
    outputs_dir = Path("outputs")
    outputs_dir.mkdir(parents=True, exist_ok=True)
    
    md_path = outputs_dir / "data_citations_table.md"
    csv_path = outputs_dir / "data_citations_table.csv"
    
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md_content)
        
    df_unique.to_csv(csv_path, index=False)
    
    print(f"✓ Saved publication quality table to: {md_path}")
    print(f"✓ Saved raw CSV table to: {csv_path}")

if __name__ == "__main__":
    main()
