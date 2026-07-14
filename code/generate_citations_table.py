#!/usr/bin/env python3
"""
generate_citations_table.py
===========================
Extracts camera trap citation data and cross-references with the processed outputs
(camera_traps_joint_metrics.csv, camera_traps_cluster_level_metrics_robust.csv,
and camera_traps_joint_detections.csv) to only include projects actually used in the
final models. Computes sample sizes: deployments provided vs. used, sorting the final
publication-quality table by used deployments.
"""

import os
import pandas as pd
from pathlib import Path

def make_markdown_table(df):
    """Generates a markdown table string from a pandas DataFrame without iterrows."""
    headers = list(df.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |"
    ]
    for row in df.to_dict("records"):
        val_str = [
            str(row[h]).replace("|", "\\|").replace("\n", " ").replace("\r", " ").strip()
            for h in headers
        ]
        lines.append("| " + " | ".join(val_str) + " |")
    return "\n".join(lines)

def main():
    print("=== Cross-referencing Camera Trap Citations with Processing Outputs ===")
    
    # 1. Load robust cluster IDs (used in the final analysis)
    robust_metrics_path = Path("outputs/camera_traps_cluster_level_metrics_robust.csv")
    if not robust_metrics_path.exists():
        print(f"Error: {robust_metrics_path} is missing. Please run the integration tests first.")
        return
    robust_clusters_df = pd.read_csv(robust_metrics_path)
    robust_cluster_ids = set(robust_clusters_df["cluster_id"].astype(str))
    print(f"Loaded {len(robust_cluster_ids)} robust cluster IDs used in the models.")

    # 2. Load joint detections to map deployments to cluster_id
    detections_path = Path("outputs/camera_traps_joint_detections.csv")
    if not detections_path.exists():
        print(f"Error: {detections_path} is missing.")
        return
    
    print("Loading event-level joint detections...")
    # Load columns needed for mapping
    det_cols = ["project_id", "deployment_id", "cluster_id"]
    det_df = pd.read_csv(detections_path, usecols=det_cols, keep_default_na=False)
    
    # Standardize column types
    det_df["project_id"] = det_df["project_id"].astype(str).str.strip()
    det_df["deployment_id"] = det_df["deployment_id"].astype(str).str.strip()
    det_df["cluster_id"] = det_df["cluster_id"].astype(str).str.strip()
    
    # Map deployments to their cluster ID
    dep_to_cluster = {}
    for dep_id, cl_id in zip(det_df["deployment_id"], det_df["cluster_id"]):
        if dep_id and cl_id:
            dep_to_cluster[dep_id] = cl_id

    # 3. Find and scan all raw packages to get provided deployments and citations
    data_dir = Path("data")
    deployments_files = list(data_dir.glob("**/deployments.csv"))
    print(f"Found {len(deployments_files)} raw deployments.csv files in data/ directory.")
    
    project_metadata = {}
    
    for dep_path in deployments_files:
        package_dir = dep_path.parent
        proj_csv_path = package_dir / "projects.csv"
        if not proj_csv_path.exists():
            continue
            
        # Parse projects.csv for name and citation
        try:
            proj_df = pd.read_csv(proj_csv_path, keep_default_na=False)
            required = ["project_id", "project_name", "data_citation"]
            if not all(col in proj_df.columns for col in required):
                col_map = {col.lower(): col for col in proj_df.columns}
                if all(req.lower() in col_map for req in required):
                    proj_df = proj_df.rename(columns={col_map[req.lower()]: req for req in required})
                else:
                    continue
            
            for row in proj_df.to_dict("records"):
                pid = str(row["project_id"]).strip()
                pname = str(row["project_name"]).strip()
                cit = str(row["data_citation"]).strip()
                
                # Determine basin based on path
                path_str = str(package_dir)
                if "Amazon" in path_str:
                    basin = "Amazon"
                elif "Congo" in path_str:
                    basin = "Congo"
                elif "SEAsia" in path_str:
                    basin = "SE Asia"
                else:
                    basin = "Other"
                    
                if pid and cit:
                    if pid not in project_metadata:
                        project_metadata[pid] = {
                            "name": pname,
                            "citation": cit,
                            "basin": basin,
                            "provided_deps": set(),
                            "used_deps": set()
                        }
        except Exception as e:
            print(f"Error reading projects.csv at {proj_csv_path}: {e}")

        # Parse deployments.csv to count provided deployments
        try:
            dep_df = pd.read_csv(dep_path, keep_default_na=False)
            # Group deployments by project_id to avoid row-by-row iterrows loops
            dep_groups = dep_df.groupby("project_id")["deployment_id"].apply(set).to_dict()
            
            for pid, deps in dep_groups.items():
                pid_str = str(pid).strip()
                if pid_str not in project_metadata:
                    continue
                    
                meta = project_metadata[pid_str]
                if "provided_deps" not in meta:
                    meta["provided_deps"] = set()
                    meta["used_deps"] = set()
                    
                cleaned_deps = {str(d).strip() for d in deps}
                meta["provided_deps"].update(cleaned_deps)
                
                # Filter used deployments using vectorized/set intersection
                used_in_file = {
                    d for d in cleaned_deps 
                    if d in dep_to_cluster and dep_to_cluster[d] in robust_cluster_ids
                }
                meta["used_deps"].update(used_in_file)
                        
        except Exception as e:
            print(f"Error reading deployments.csv at {dep_path}: {e}")

    # 4. Filter and compile the final citations list
    records = []
    total_provided_deps = 0
    total_used_deps = 0
    
    for pid, meta in project_metadata.items():
        # Only include if deployments were actually used in final robust models
        used_deps_set = meta.get("used_deps", set())
        used_deps_count = len(used_deps_set)
        
        if used_deps_count > 0:
            provided_deps_count = len(meta.get("provided_deps", set()))
            
            total_provided_deps += provided_deps_count
            total_used_deps += used_deps_count
            
            records.append({
                "Basin": meta["basin"],
                "Project Name": meta["name"],
                "Project ID": pid,
                "Deployments (Provided)": provided_deps_count,
                "Deployments (Used)": used_deps_count,
                "Citation": meta["citation"]
            })
            
    if not records:
        print("No active project citations found! Exiting.")
        return
        
    df_citations = pd.DataFrame(records)
    
    # 5. Sort by Deployments Used (descending), then Project Name
    df_citations = df_citations.sort_values(
        by=["Deployments (Used)", "Project Name"],
        ascending=[False, True]
    ).reset_index(drop=True)
    
    print(f"\nFinal statistics across used datasets:")
    print(f"  Unique Projects Used: {len(df_citations)} (out of {len(project_metadata)} total provided)")
    print(f"  Deployments Provided: {total_provided_deps} | Used: {total_used_deps}")
    
    # 6. Generate Markdown publication-quality document
    md_content = "# Camera Trap Data Citations & Sample Sizes\n\n"
    md_content += "This table lists all camera trap datasets actually used in the standing mammal biomass calibration models, "
    md_content += "sorted by the number of used deployments. It compares the number of deployments provided by each project "
    md_content += "to those retained in robust clusters for final RAI calibrations.\n\n"
    
    # Create a nice summary header block
    md_content += "### Data Summary\n"
    md_content += f"- **Total Projects Used**: {len(df_citations)}\n"
    md_content += f"- **Deployments Provided**: {total_provided_deps:,} | **Deployments Used**: {total_used_deps:,}\n\n"
    
    # Append the main markdown table
    md_content += make_markdown_table(df_citations)
    md_content += "\n"
    
    # Write outputs
    outputs_dir = Path("outputs")
    outputs_dir.mkdir(parents=True, exist_ok=True)
    
    md_path = outputs_dir / "data_citations_table.md"
    csv_path = outputs_dir / "data_citations_table.csv"
    
    with open(md_path, "w", encoding="utf-8") as f:
        f.write(md_content)
        
    df_citations.to_csv(csv_path, index=False)
    
    print(f"\n✓ Saved publication quality table to: {md_path}")
    print(f"✓ Saved raw CSV table to: {csv_path}")

if __name__ == "__main__":
    main()
