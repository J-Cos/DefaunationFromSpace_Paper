#!/usr/bin/env python3
"""
process_camera_traps.py
=======================
Ingests both Amazon Basin and Congo Basin Wildlife Insights datasets,
collapses image-level data to independent events (30-min threshold),
matches taxa to EltonTraits body-mass databases, calculates biophysical
metrics (B_H and M_H indices) at the deployment level, and exports
joint tidy datasets.

Input:  - data/CongoCameraTrapping/wildlife-insights_* (Congo Basin packages)
        - data/AmazonCameraTrapping/wildlife-insights_* (Amazon Basin packages)
        - data/trait_databases/MamFuncDat.txt           (EltonTraits Mammals)
        - data/trait_databases/BirdFuncDat.txt          (EltonTraits Birds)

Output: - outputs/camera_traps_joint_detections.csv     (Event-level joint detections)
        - outputs/camera_traps_joint_metrics.csv        (Deployment-level joint metrics)
"""

import os
import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from datetime import datetime

# ── Paths ───────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

CONGO_DIR = DATA_DIR / "CongoCameraTrapping"
AMAZON_DIR = DATA_DIR / "AmazonCameraTrapping"
TRAIT_DIR = DATA_DIR / "trait_databases"

# ── Allometric & Biophysical Constants ──────────────────────────────────────
# Day range (km/day) from body mass (kg): Carbone et al. 2005
DAY_RANGE_COEFF = 1.2
DAY_RANGE_EXP = 0.26

# Field metabolic rate (Watts) from body mass (kg): Nagy 2005
FMR_COEFF = 10.0  # Watts
FMR_EXP = 0.75

# ── Fallback Body Masses (kg) by Family for Taxa not in EltonTraits ─────────
# Updated to cover both African and South American tropical forest species
FALLBACK_FAMILY_MASS_KG = {
    # African / Global Fallbacks
    "Viverridae": 2.5,
    "Herpestidae": 1.5,
    "Muridae": 0.04,
    "Sciuridae": 0.3,
    "Nesomyidae": 0.8,
    "Columbidae": 0.2,
    "Numididae": 1.3,
    "Phasianidae": 1.0,
    "Pittidae": 0.08,
    "Megapodiidae": 1.5,
    "Cracidae": 2.0,
    "Cervidae": 30.0,
    "Tragulidae": 3.0,
    "Canidae": 10.0,
    "Felidae": 15.0,
    "Mustelidae": 2.0,
    "Procyonidae": 5.0,
    "Manidae": 8.0,
    "Didelphidae": 0.5,
    "Hystricidae": 3.0,
    "Elephantidae": 4000.0,
    "Suidae": 60.0,
    "Bovidae": 20.0,
    "Cercopithecidae": 7.0,
    "Hominidae": 70.0,
    
    # Neotropical Fallbacks (Amazon)
    "Myrmecophagidae": 15.0,   # Sloths and anteaters
    "Bradypodidae": 4.5,
    "Megalonychidae": 6.0,
    "Dasypodidae": 4.5,       # Armadillos
    "Chlamyphoridae": 4.0,
    "Tayassuidae": 25.0,      # Peccaries
    "Cuniculidae": 8.0,       # Pacas
    "Dasyproctidae": 3.5,     # Agoutis
    "Echimyidae": 0.4,        # Spiny rats
    "Cebidae": 3.0,           # New World monkeys
    "Atelidae": 8.0,
    "Callitrichidae": 0.4,
    "Tapiridae": 200.0,       # Tapirs
    "Tinamidae": 1.0,         # Tinamous
    "Psophiidae": 1.0,        # Trumpeters
    "Cariamidae": 1.5,        # Seriemas
    "Ramphastidae": 0.6,      # Toucans
}

DOMESTIC_SPECIES = {
    ("Bos", "taurus"),       # Domestic Cattle
    ("Canis", "familiaris"), # Domestic Dog
    ("Felis", "catus"),      # Domestic Cat
    ("Sus", "domesticus"),   # Domestic Pig
    ("Capra", "hircus"),     # Domestic Goat
    ("Ovis", "aries"),       # Domestic Sheep
    ("Equus", "caballus"),   # Domestic Horse
}

EXCLUDE_COMMON_NAMES = {
    "Blank", "No CV Result", "Human", "Homo Species",
    "Setup", "Calibration", "Unknown", "Vehicle",
}


# ── Ingestion Helpers ────────────────────────────────────────────────────────

def discover_wi_packages(region_dir: Path) -> list[Path]:
    """Find all unzipped WI data package directories under a folder."""
    packages = []
    if not region_dir.exists():
        return packages
    for d in sorted(region_dir.iterdir()):
        if d.is_dir() and (d / "deployments.csv").exists():
            packages.append(d)
    return packages


def load_deployments(package_dir: Path) -> pd.DataFrame:
    """Load and clean the deployments table from a single WI package."""
    df = pd.read_csv(package_dir / "deployments.csv")
    # Clean column names
    df.columns = df.columns.str.strip().str.lower()
    
    # Parse dates
    df["start_date"] = pd.to_datetime(df["start_date"], errors="coerce")
    df["end_date"] = pd.to_datetime(df["end_date"], errors="coerce")
    
    # Compute trap-days
    df["trap_days"] = (df["end_date"] - df["start_date"]).dt.total_seconds() / 86400
    
    # Drop deployments with missing or invalid dates or locations
    df = df.dropna(subset=["start_date", "end_date", "longitude", "latitude"])
    df = df[df["trap_days"] > 0].copy()
    df = df.drop_duplicates(subset=["project_id", "deployment_id"]).copy()
    
    cols = ["project_id", "deployment_id", "longitude", "latitude",
            "start_date", "end_date", "trap_days"]
    # Check if optional project name/subproject exists
    for col in ["project_name", "subproject_name", "feature_type"]:
        if col in df.columns:
            cols.append(col)
    
    return df[cols]


def _parse_wi_timestamp(ts_series: pd.Series) -> pd.Series:
    """Parse Wildlife Insights timestamps supporting standard ISO or GMT string formats."""
    result = pd.to_datetime(ts_series, errors="coerce", utc=True)
    mask = result.isna() & ts_series.notna()
    if mask.any():
        cleaned = ts_series[mask].astype(str).str.replace(
            r"\s*\(.*\)\s*$", "", regex=True
        )
        result[mask] = pd.to_datetime(cleaned, errors="coerce", utc=True)
    return result


def _collapse_to_independent_events(
    df: pd.DataFrame,
    threshold_minutes: float = 30.0,
) -> pd.DataFrame:
    """
    Collapse image-level detections into sequence-level independent events.
    Groups images of the same taxon at the same deployment within a threshold.
    """
    keep_cols = ["project_id", "deployment_id", "genus", "species",
                 "common_name", "class", "order", "family",
                 "number_of_objects", "timestamp"]
    keep_cols = [c for c in keep_cols if c in df.columns]

    # Split by sequence_id presence
    has_seq_id = "sequence_id" in df.columns
    if has_seq_id:
        with_seq = df[df["sequence_id"].notna()].copy()
        without_seq = df[df["sequence_id"].isna()].copy()
    else:
        with_seq = pd.DataFrame(columns=df.columns)
        without_seq = df.copy()

    events = []

    if len(with_seq) > 0:
        # Collapse using sequence_id, grouping by taxon to preserve separate species/families inside the sequence
        grp_cols_seq = ["project_id", "deployment_id", "sequence_id"]
        for c in ["class", "order", "family", "genus", "species"]:
            if c in with_seq.columns:
                grp_cols_seq.append(c)

        seq_events = (with_seq
                      .groupby(grp_cols_seq, dropna=False)
                      .agg(
                          common_name=("common_name", "first"),
                          number_of_objects=("number_of_objects", "max"),
                          timestamp=("timestamp", "first"),
                      )
                      .reset_index()
                      .drop(columns=["sequence_id"]))
        events.append(seq_events)

    if len(without_seq) > 0 and "timestamp" in without_seq.columns:
        # Collapse using temporal threshold
        without_seq["_ts"] = _parse_wi_timestamp(without_seq["timestamp"])
        
        # Sort and group by deployment x taxon hierarchy to keep different families/genera separate
        grp_cols = ["project_id", "deployment_id"]
        for c in ["class", "order", "family", "genus", "species"]:
            if c in without_seq.columns:
                grp_cols.append(c)
                
        without_seq = without_seq.sort_values(grp_cols + ["_ts"])
        
        without_seq["_gap"] = (
            without_seq.groupby(grp_cols, dropna=False)["_ts"]
            .diff()
            .dt.total_seconds()
            .fillna(threshold_minutes * 60 + 1)
        )
        without_seq["_new_event"] = without_seq["_gap"] > (threshold_minutes * 60)
        without_seq["_event_id"] = without_seq.groupby(grp_cols, dropna=False)["_new_event"].cumsum()

        agg_dict = {
            "common_name": ("common_name", "first"),
            "number_of_objects": ("number_of_objects", "max"),
            "timestamp": ("timestamp", "first"),
        }

        temporal_events = (without_seq
                           .groupby(grp_cols + ["_event_id"], dropna=False)
                           .agg(**agg_dict)
                           .reset_index()
                           .drop(columns=["_event_id"]))
        events.append(temporal_events)
    elif len(without_seq) > 0:
        events.append(without_seq[keep_cols])

    if not events:
        return pd.DataFrame(columns=keep_cols)

    result = pd.concat(events, ignore_index=True)
    out_cols = [c for c in keep_cols if c in result.columns]
    return result[out_cols]


def load_images(package_dir: Path, independence_threshold_min: float = 30.0) -> pd.DataFrame:
    """Load and collapse image-level or sequence-level records from a WI package."""
    img_path = package_dir / "images.csv"
    seq_path = package_dir / "sequences.csv"

    if seq_path.exists():
        # Sequences exists: pre-collapsed event-level data
        df = pd.read_csv(seq_path)
        df.columns = df.columns.str.strip().str.lower()

        if "group_size" in df.columns and "number_of_objects" not in df.columns:
            df["number_of_objects"] = pd.to_numeric(
                df["group_size"], errors="coerce"
            ).fillna(1).astype(int)

        if "timestamp" not in df.columns and "start_time" in df.columns:
            df["timestamp"] = df["start_time"]

        keep = ["project_id", "deployment_id", "genus", "species",
                "common_name", "class", "order", "family",
                "number_of_objects", "timestamp"]
        keep = [c for c in keep if c in df.columns]
        df = df[keep].copy()

        df["number_of_objects"] = pd.to_numeric(
            df.get("number_of_objects", 1), errors="coerce"
        ).fillna(1).astype(int)

        return df

    elif img_path.exists():
        # Images exists: requires event collapsing
        df = pd.read_csv(img_path)
        df.columns = df.columns.str.strip().str.lower()

        df["number_of_objects"] = pd.to_numeric(
            df.get("number_of_objects", 1), errors="coerce"
        ).fillna(1).astype(int)

        n_before = len(df)
        df = _collapse_to_independent_events(df, independence_threshold_min)
        n_after = len(df)
        if n_before != n_after:
            print(f"    Collapsed {n_before} images -> {n_after} independent events "
                  f"({n_before/max(n_after,1):.1f}x reduction)")

        return df
    else:
        raise FileNotFoundError(f"No images.csv or sequences.csv in {package_dir}")


def safe_str(val) -> str:
    """Helper to cleanly stringify dataframe elements, replacing float NaN with empty string."""
    if pd.isna(val):
        return ""
    return str(val).strip()


def assign_taxon_quality(row: pd.Series) -> str:
    """Classify taxonomic resolution quality and flag domestic/human/blanks."""
    cn = safe_str(row.get("common_name", ""))
    genus = safe_str(row.get("genus", ""))
    species = safe_str(row.get("species", ""))
    family = safe_str(row.get("family", ""))

    cn_lower = cn.lower()
    genus_lower = genus.lower()
    family_lower = family.lower()

    # Check for human detections
    is_human = (
        genus_lower == "homo" or
        "human" in cn_lower or
        "trapper" in cn_lower or
        "researcher" in cn_lower or
        "hunter" in cn_lower or
        "pedestrian" in cn_lower or
        "rider" in cn_lower or
        "biker" in cn_lower or
        "maintenance" in cn_lower
    )
    if is_human:
        return "human"

    if cn in EXCLUDE_COMMON_NAMES or cn == "" or genus_lower in ["setup", "unknown", "vehicle"]:
        return "blank"

    if (genus, species) in DOMESTIC_SPECIES:
        return "domestic"

    if genus and species and genus != "" and species != "" and species.lower() not in ("sp", "sp."):
        return "species"
    if genus and genus != "" and genus_lower != "unknown":
        return "genus"
    if family and family != "" and family_lower != "unknown":
        return "family"
    return "higher"



def compute_detection_rates(deployments: pd.DataFrame,
                             images: pd.DataFrame) -> pd.DataFrame:
    """Group by taxon, count detections, and calculate RAI per 100 trap-days."""
    images["taxon_quality"] = images.apply(assign_taxon_quality, axis=1)

    images["taxon_key"] = images.apply(
        lambda r: f"{safe_str(r.get('genus'))}_{safe_str(r.get('species'))}"
        if r["taxon_quality"] in ("species", "genus")
        else f"{safe_str(r.get('family'))}_{safe_str(r.get('common_name'))}",
        axis=1
    )

    det = (images
           .groupby(["project_id", "deployment_id", "taxon_key",
                     "class", "order", "family", "genus", "species",
                     "taxon_quality"], dropna=False)
           .agg(n_detections=("number_of_objects", "sum"),
                common_name=("common_name", "first"))
           .reset_index())

    det = det.merge(
        deployments,
        on=["project_id", "deployment_id"],
        how="right"
    )

    det["n_detections"] = det["n_detections"].fillna(0)
    det["taxon_key"] = det["taxon_key"].fillna("None_None")
    det["taxon_quality"] = det["taxon_quality"].fillna("blank")
    det["RAI"] = (det["n_detections"] / det["trap_days"]) * 100

    return det.sort_values(["project_id", "deployment_id", "taxon_key"])


# ── Trait Ingestion & Allometric Metrics ─────────────────────────────────────

def load_eltontraits(trait_dir: Path) -> pd.DataFrame:
    """Load EltonTraits mammal and bird body-mass databases."""
    frames = []

    # Mammals
    mam_path = trait_dir / "MamFuncDat.txt"
    if mam_path.exists():
        mam = pd.read_csv(mam_path, sep="\t", encoding="latin-1")
        if "Scientific" in mam.columns:
            mam[["genus", "species"]] = mam["Scientific"].str.strip().str.split(" ", n=1, expand=True)
            mass_col = [c for c in mam.columns if "BodyMass" in c and "Value" in c]
            if mass_col:
                mam["body_mass_g"] = pd.to_numeric(mam[mass_col[0]], errors="coerce")
                mam["body_mass_kg"] = mam["body_mass_g"] / 1000
                mam["class"] = "Mammalia"
                mam["taxon_source"] = "EltonTraits_mammal"
                frames.append(mam[["genus", "species", "body_mass_kg", "class", "taxon_source"]])
                print(f"  Loaded {len(mam)} mammal species from EltonTraits")

    # Birds
    bird_path = trait_dir / "BirdFuncDat.txt"
    if bird_path.exists():
        bird = pd.read_csv(bird_path, sep="\t", encoding="latin-1")
        if "Scientific" in bird.columns:
            bird[["genus", "species"]] = bird["Scientific"].str.strip().str.split(" ", n=1, expand=True)
            mass_col = [c for c in bird.columns if "BodyMass" in c and "Value" in c]
            if mass_col:
                bird["body_mass_g"] = pd.to_numeric(bird[mass_col[0]], errors="coerce")
                bird["body_mass_kg"] = bird["body_mass_g"] / 1000
                bird["class"] = "Aves"
                bird["taxon_source"] = "EltonTraits_bird"
                frames.append(bird[["genus", "species", "body_mass_kg", "class", "taxon_source"]])
                print(f"  Loaded {len(bird)} bird species from EltonTraits")

    if not frames:
        print("  WARNING: No EltonTraits files found. Using fallback masses only.")
        return pd.DataFrame(columns=["genus", "species", "body_mass_kg", "class", "taxon_source"])

    return pd.concat(frames, ignore_index=True).dropna(subset=["body_mass_kg"])


def match_body_mass(det: pd.DataFrame, traits: pd.DataFrame) -> pd.DataFrame:
    """Matches body mass at exact species binomial, genus median, or family fallback levels."""
    det = det.copy()
    det["body_mass_kg"] = np.nan
    det["mass_match_level"] = "unmatched"

    if len(traits) > 0:
        species_mass = traits.groupby(["genus", "species"])["body_mass_kg"].median().to_dict()
        genus_mass = traits.groupby("genus")["body_mass_kg"].median().to_dict()
    else:
        species_mass = {}
        genus_mass = {}

    for idx, row in det.iterrows():
        g = safe_str(row.get("genus", ""))
        s = safe_str(row.get("species", ""))
        fam = safe_str(row.get("family", ""))

        # 1. Exact species match
        if g and s and (g, s) in species_mass:
            det.at[idx, "body_mass_kg"] = species_mass[(g, s)]
            det.at[idx, "mass_match_level"] = "species"
            continue

        # 2. Genus median fallback
        if g and g in genus_mass:
            det.at[idx, "body_mass_kg"] = genus_mass[g]
            det.at[idx, "mass_match_level"] = "genus"
            continue

        # 3. Family fallback
        if fam and fam in FALLBACK_FAMILY_MASS_KG:
            det.at[idx, "body_mass_kg"] = FALLBACK_FAMILY_MASS_KG[fam]
            det.at[idx, "mass_match_level"] = "family_fallback"
            continue

    return det


def compute_site_heterotroph_metrics(det: pd.DataFrame) -> pd.DataFrame:
    """Calculate site-level heterotroph biomass and metabolism indices using Carbone & Nagy corrections."""
    usable = det[
        (det["taxon_quality"].isin(["species", "genus", "family"])) &
        (~det["taxon_quality"].isin(["blank", "human", "domestic"])) &
        (det["body_mass_kg"].notna()) &
        (det["body_mass_kg"] > 0)
    ].copy()

    # Day-range correction
    usable["day_range_km"] = DAY_RANGE_COEFF * (usable["body_mass_kg"] ** DAY_RANGE_EXP)
    usable["corrected_RAI"] = usable["RAI"] / usable["day_range_km"]

    # Biophysical contributions
    usable["biomass_contrib"] = usable["corrected_RAI"] * usable["body_mass_kg"]
    usable["metabolism_contrib"] = usable["corrected_RAI"] * FMR_COEFF * (usable["body_mass_kg"] ** FMR_EXP)

    # Large mammal contributions
    usable["biomass_contrib_gt50"] = np.where(usable["body_mass_kg"] > 50.0, usable["biomass_contrib"], 0.0)
    usable["biomass_contrib_gt100"] = np.where(usable["body_mass_kg"] > 100.0, usable["biomass_contrib"], 0.0)

    agg_dict = {
        "n_species": ("taxon_key", "nunique"),
        "n_detections_total": ("n_detections", "sum"),
        "B_H_index": ("biomass_contrib", "sum"),
        "M_H_index": ("metabolism_contrib", "sum"),
        "B_H_gt50": ("biomass_contrib_gt50", "sum"),
        "B_H_gt100": ("biomass_contrib_gt100", "sum"),
        "dominant_species": ("common_name", lambda x: str(x.dropna().mode().iloc[0]) if len(x.dropna().mode()) > 0 else ""),
        "species_list": ("common_name", lambda x: "; ".join(sorted(x.dropna().astype(str).unique()))),
    }

    site = (usable
            .groupby(["project_id", "deployment_id"])
            .agg(**agg_dict)
            .reset_index())

    # Calculate fraction
    site["megafauna_fraction"] = np.where(site["B_H_index"] > 0, (site["B_H_gt50"] / site["B_H_index"]) * 100, 0.0)

    # Retain all deployments (even those with zero wild species)
    meta_cols = ["project_id", "deployment_id", "longitude", "latitude", "trap_days"]
    if "project_name" in det.columns:
        meta_cols.append("project_name")
    
    all_deps = det[meta_cols].drop_duplicates()
    site = all_deps.merge(site, on=["project_id", "deployment_id"], how="left")

    site["n_species"] = site["n_species"].fillna(0).astype(int)
    site["n_detections_total"] = site["n_detections_total"].fillna(0).astype(int)
    site["B_H_index"] = site["B_H_index"].fillna(0.0)
    site["M_H_index"] = site["M_H_index"].fillna(0.0)
    site["B_H_gt50"] = site["B_H_gt50"].fillna(0.0)
    site["B_H_gt100"] = site["B_H_gt100"].fillna(0.0)
    site["megafauna_fraction"] = site["megafauna_fraction"].fillna(0.0)
    site["dominant_species"] = site["dominant_species"].fillna("")
    site["species_list"] = site["species_list"].fillna("")

    return site


# ── Pipeline Driver ──────────────────────────────────────────────────────────

def process_region(region_dir: Path, region_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Discover, load, and compute metrics for all packages in a given region."""
    packages = discover_wi_packages(region_dir)
    print(f"\nProcessing {region_name} Basin ({len(packages)} packages):")
    
    all_deployments = []
    all_events = []
    
    for pkg in packages:
        print(f"  Ingesting {pkg.name}...")
        dep = load_deployments(pkg)
        evt = load_images(pkg, independence_threshold_min=30.0)
        
        # Merge project name from projects.csv if available
        proj_csv = pkg / "projects.csv"
        if proj_csv.exists():
            try:
                proj_df = pd.read_csv(proj_csv)
                proj_df.columns = proj_df.columns.str.strip().str.lower()
                if "project_name" in proj_df.columns and "project_id" in proj_df.columns:
                    name_map = proj_df.set_index("project_id")["project_name"].to_dict()
                    dep["project_name"] = dep["project_id"].map(name_map)
            except Exception as e:
                print(f"    Warning: Could not parse projects.csv: {e}")
        
        if "project_name" not in dep.columns:
            dep["project_name"] = pkg.name
            
        all_deployments.append(dep)
        all_events.append(evt)
        print(f"    - {len(dep)} deployments, {len(evt)} independent events")

    if not all_deployments:
        print(f"  No packages found for {region_name} Basin.")
        return pd.DataFrame(), pd.DataFrame()

    deployments = pd.concat(all_deployments, ignore_index=True)
    deployments = deployments.drop_duplicates(subset=["project_id", "deployment_id"])

    events = pd.concat(all_events, ignore_index=True)
    dedup_cols = ["project_id", "deployment_id", "genus", "species", "timestamp", "number_of_objects"]
    dedup_cols = [c for c in dedup_cols if c in events.columns]
    events = events.drop_duplicates(subset=dedup_cols)

    print(f"  Total {region_name} clean deployments: {len(deployments)}")
    print(f"  Total {region_name} clean events:      {len(events)}")

    det = compute_detection_rates(deployments, events)
    det["region"] = region_name
    return det, events


def main():
    parser = argparse.ArgumentParser(
        description="Process joint Congo and Amazon camera trap datasets with EltonTraits allometric scaling."
    )
    parser.add_argument("--congo-dir", type=Path, default=CONGO_DIR)
    parser.add_argument("--amazon-dir", type=Path, default=AMAZON_DIR)
    parser.add_argument("--trait-dir", type=Path, default=TRAIT_DIR)
    parser.add_argument("--out-detections", type=Path, default=OUTPUT_DIR / "camera_traps_joint_detections.csv")
    parser.add_argument("--out-metrics", type=Path, default=OUTPUT_DIR / "camera_traps_joint_metrics.csv")
    args = parser.parse_args()

    # Load traits database
    print("Loading EltonTraits body-mass database...")
    traits = load_eltontraits(args.trait_dir)

    # Process Congo
    congo_det, congo_raw_events = process_region(args.congo_dir, "Congo")

    # Process Amazon
    amazon_det, amazon_raw_events = process_region(args.amazon_dir, "Amazon")

    # Combine detections
    print("\nMerging datasets...")
    joint_det = pd.concat([congo_det, amazon_det], ignore_index=True)

    # Match body masses
    print("Matching body masses...")
    joint_det = match_body_mass(joint_det, traits)
    
    match_summary = joint_det["mass_match_level"].value_counts()
    print("Match quality:")
    for lvl, cnt in match_summary.items():
        print(f"  {lvl:18s}: {cnt:5d} rows")

    # Save event-level detections
    joint_det.to_csv(args.out_detections, index=False)
    print(f"Saved joint detections (event-level) to {args.out_detections}")

    # Compute site-level metrics
    print("\nComputing site-level heterotroph biomass and metabolism indices...")
    joint_metrics = compute_site_heterotroph_metrics(joint_det)
    joint_metrics["region"] = joint_metrics["project_id"].map(
        joint_det.drop_duplicates("project_id").set_index("project_id")["region"]
    )
    
    # Save deployment-level metrics
    joint_metrics.to_csv(args.out_metrics, index=False)
    print(f"Saved joint metrics (deployment-level) to {args.out_metrics}")

    # Regional summary statistics
    print(f"\n{'='*60}")
    print("SUMMARY")
    print(f"{'='*60}")
    for region in ["Congo", "Amazon"]:
        sub = joint_metrics[joint_metrics["region"] == region]
        print(f"{region} Basin:")
        print(f"  Number of deployments: {len(sub)}")
        print(f"  Unique projects:        {sub['project_name'].nunique()}")
        print(f"  Total species detected: {sub['n_species'].max()}")
        print(f"  Mean species/site:      {sub['n_species'].mean():.1f}")
        print(f"  Biomass Index range:    {sub['B_H_index'].min():.2f} – {sub['B_H_index'].max():.2f}")
        print(f"  Metabolism Index range: {sub['M_H_index'].min():.2f} – {sub['M_H_index'].max():.2f}")
        print(f"  Median Biomass >50kg:   {sub['B_H_gt50'].median():.2f}")
        print(f"  Median Biomass >100kg:  {sub['B_H_gt100'].median():.2f}")
        print(f"  Mean Megafauna %:       {sub['megafauna_fraction'].mean():.1f}%")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
