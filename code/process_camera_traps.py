#!/usr/bin/env python3
"""
process_camera_traps.py
=======================
Ingests Congo, Amazon, and Southeast Asia camera trapping data,
collapses image-level data to independent events (30-min threshold),
matches taxa to EltonTraits body-mass databases, performs 11.1km
spatial clustering, filters clusters by GEDI grid overlap, computes
cluster-level Relative Abundance Indices (RAI), Biomass Index (B_H),
Metabolism Index (M_H), temporal weights, and taxonomic keep proportions.
Saves all CSV and GeoJSON spatial products.

Input:  - data/CongoCameraTrapping/wildlife-insights_* (Congo Basin packages)
        - data/AmazonCameraTrapping/wildlife-insights_* (Amazon Basin packages)
        - data/SEAsiaCameraTrapping/wildlife-insights_* (SE Asia packages)
        - data/trait_databases/MamFuncDat.txt           (EltonTraits Mammals)
        - data/trait_databases/BirdFuncDat.txt          (EltonTraits Birds)

Output: - outputs/camera_traps_joint_detections.csv     (Event-level joint detections with cluster assignments)
        - outputs/camera_traps_joint_metrics.csv        (Deployment-level joint metrics)
        - outputs/camera_traps_cluster_level_metrics.csv (GEDI-valid cluster-level metrics)
        - outputs/camera_traps_cluster_level_metrics_robust.csv (Robust cluster-level metrics, trap_days >= min_trap_days from config)
        - outputs/camera_traps_robust_detections.csv    (Event-level detections in robust clusters only)
        - outputs/camera_traps_robust_buffered_mcps.geojson (Robust cluster buffered MCP geometries with attributes)
"""

import os
import json
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from datetime import datetime

# Optional dependencies for spatial operations (will be checked in main)
try:
    from scipy.spatial.distance import pdist
    from scipy.cluster.hierarchy import linkage, fcluster
    from shapely.geometry import MultiPoint, mapping
    from shapely.ops import transform
    from pyproj import Transformer
    import rasterio
    from rasterio.mask import mask
    SPATIAL_LIBS_AVAILABLE = True
except ImportError:
    SPATIAL_LIBS_AVAILABLE = False

# ── Paths ───────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
DATA_DIR = PROJECT_DIR / "data"
OUTPUT_DIR = PROJECT_DIR / "outputs"
OUTPUT_DIR.mkdir(exist_ok=True)

CONGO_DIR = DATA_DIR / "CongoCameraTrapping"
AMAZON_DIR = DATA_DIR / "AmazonCameraTrapping"
SEASIA_DIR = DATA_DIR / "SEAsiaCameraTrapping"
TRAIT_DIR = DATA_DIR / "trait_databases"

# ── Load Config File ────────────────────────────────────────────────────────
CONFIG_PATH = SCRIPT_DIR / "config.json"
if CONFIG_PATH.exists():
    with open(CONFIG_PATH, "r") as f:
        CONFIG = json.load(f)
else:
    # Fallback default values
    CONFIG = {
        "allometric_scaling": {
            "day_range_coeff": 1.2,
            "day_range_exp": 0.26,
            "fmr_coeff": 5.7,
            "fmr_exp": 0.75
        },
        "clustering": {
            "threshold_km": 11.1,
            "buffer_meters": 5500,
            "min_trap_days": 10
        },
        "temporal_decay": {
            "brackets": [1.0, 6.0, 11.0],
            "weights": [1.0, 0.5, 0.25, 0.1]
        }
    }

DAY_RANGE_COEFF = CONFIG["allometric_scaling"]["day_range_coeff"]
DAY_RANGE_EXP = CONFIG["allometric_scaling"]["day_range_exp"]
FMR_COEFF = CONFIG["allometric_scaling"]["fmr_coeff"]
FMR_EXP = CONFIG["allometric_scaling"]["fmr_exp"]
CLUSTER_THRESHOLD_KM = CONFIG["clustering"]["threshold_km"]
BUFFER_METERS = CONFIG["clustering"]["buffer_meters"]
MIN_TRAP_DAYS = CONFIG["clustering"]["min_trap_days"]

# ── Fallback Body Masses (kg) by Family for Taxa not in EltonTraits ─────────
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
    "Myrmecophagidae": 15.0,   # Anteaters
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
    for path in sorted(region_dir.rglob("deployments.csv")):
        if path.is_file():
            packages.append(path.parent)
    return packages


def load_deployments(package_dir: Path) -> pd.DataFrame:
    """Load and clean the deployments table from a single WI package."""
    df = pd.read_csv(package_dir / "deployments.csv")
    df.columns = df.columns.str.strip().str.lower()
    
    df["start_date"] = pd.to_datetime(df["start_date"], errors="coerce")
    df["end_date"] = pd.to_datetime(df["end_date"], errors="coerce")
    df["trap_days"] = (df["end_date"] - df["start_date"]).dt.total_seconds() / 86400
    
    df = df.dropna(subset=["start_date", "end_date", "longitude", "latitude"])
    df = df[df["trap_days"] > 0].copy()
    df = df.drop_duplicates(subset=["project_id", "deployment_id"]).copy()
    
    cols = ["project_id", "deployment_id", "longitude", "latitude",
            "start_date", "end_date", "trap_days"]
    for col in ["project_name", "subproject_name", "feature_type"]:
        if col in df.columns:
            cols.append(col)
    
    return df[cols]


def _parse_wi_timestamp(ts_series: pd.Series) -> pd.Series:
    """Parse Wildlife Insights timestamps supporting ISO or GMT string formats."""
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
    """Collapse image-level detections into sequence-level independent events."""
    keep_cols = ["project_id", "deployment_id", "genus", "species",
                 "common_name", "class", "order", "family",
                 "number_of_objects", "timestamp"]
    keep_cols = [c for c in keep_cols if c in df.columns]

    has_seq_id = "sequence_id" in df.columns
    if has_seq_id:
        with_seq = df[df["sequence_id"].notna()].copy()
        without_seq = df[df["sequence_id"].isna()].copy()
    else:
        with_seq = pd.DataFrame(columns=df.columns)
        without_seq = df.copy()

    events = []

    if len(with_seq) > 0:
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
        without_seq["_ts"] = _parse_wi_timestamp(without_seq["timestamp"])
        
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
    img_parts = sorted(package_dir.glob("images_*.csv"))

    if seq_path.exists():
        df = pd.read_csv(seq_path, low_memory=False)
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

    elif img_path.exists() or len(img_parts) > 0:
        if img_path.exists():
            df = pd.read_csv(img_path, low_memory=False)
        else:
            print(f"    Found split image files: {[p.name for p in img_parts]}")
            df = pd.concat([pd.read_csv(p, low_memory=False) for p in img_parts], ignore_index=True)

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
        raise FileNotFoundError(f"No images.csv, images_*.csv or sequences.csv in {package_dir}")


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
    """Matches body mass at exact species binomial, genus median, or family fallback levels.

    Uses vectorized lookups instead of row-wise iteration for performance.
    Priority: species-level > genus-level > family fallback.
    """
    det = det.copy()
    det["body_mass_kg"] = np.nan
    det["mass_match_level"] = "unmatched"

    # Ensure string columns for lookup keys
    det["_g"] = det["genus"].apply(safe_str)
    det["_s"] = det["species"].apply(safe_str)
    det["_f"] = det["family"].apply(safe_str)

    if len(traits) > 0:
        species_mass = traits.groupby(["genus", "species"])["body_mass_kg"].median()
        genus_mass = traits.groupby("genus")["body_mass_kg"].median()
    else:
        species_mass = pd.Series(dtype=float)
        genus_mass = pd.Series(dtype=float)

    # 1. Species-level match (vectorized via MultiIndex map)
    species_key = list(zip(det["_g"], det["_s"]))
    species_idx = pd.MultiIndex.from_tuples(species_key, names=["genus", "species"])
    sp_vals = species_idx.map(species_mass.to_dict().get)
    sp_mask = pd.notna(sp_vals) & (det["_g"] != "") & (det["_s"] != "")
    det.loc[sp_mask, "body_mass_kg"] = sp_vals[sp_mask]
    det.loc[sp_mask, "mass_match_level"] = "species"

    # 2. Genus-level match (only where not yet matched)
    unmatched = det["mass_match_level"] == "unmatched"
    genus_vals = det.loc[unmatched, "_g"].map(genus_mass.to_dict())
    g_mask = unmatched & genus_vals.notna() & (det["_g"] != "")
    det.loc[g_mask, "body_mass_kg"] = genus_vals[g_mask]
    det.loc[g_mask, "mass_match_level"] = "genus"

    # 3. Family fallback (only where still unmatched)
    unmatched = det["mass_match_level"] == "unmatched"
    fam_vals = det.loc[unmatched, "_f"].map(FALLBACK_FAMILY_MASS_KG)
    f_mask = unmatched & fam_vals.notna() & (det["_f"] != "")
    det.loc[f_mask, "body_mass_kg"] = fam_vals[f_mask]
    det.loc[f_mask, "mass_match_level"] = "family_fallback"

    # Clean up temporary columns
    det.drop(columns=["_g", "_s", "_f"], inplace=True)

    return det


def compute_site_heterotroph_metrics(det: pd.DataFrame) -> pd.DataFrame:
    """Calculate site-level heterotroph biomass and metabolism indices at deployment level."""
    usable = det[
        (det["taxon_quality"].isin(["species", "genus", "family"])) &
        (~det["taxon_quality"].isin(["blank", "human", "domestic"])) &
        (det["body_mass_kg"].notna()) &
        (det["body_mass_kg"] > 0)
    ].copy()

    usable["day_range_km"] = DAY_RANGE_COEFF * (usable["body_mass_kg"] ** DAY_RANGE_EXP)
    usable["corrected_RAI"] = usable["RAI"] / usable["day_range_km"]

    usable["biomass_contrib"] = usable["corrected_RAI"] * usable["body_mass_kg"]
    usable["metabolism_contrib"] = usable["corrected_RAI"] * FMR_COEFF * (usable["body_mass_kg"] ** FMR_EXP)

    usable["biomass_contrib_gt50"] = np.where(usable["body_mass_kg"] > 50.0, usable["biomass_contrib"], 0.0)
    usable["biomass_contrib_gt100"] = np.where(usable["body_mass_kg"] > 100.0, usable["biomass_contrib"], 0.0)
    usable["biomass_contrib_gt1000"] = np.where(usable["body_mass_kg"] > 1000.0, usable["biomass_contrib"], 0.0)

    agg_dict = {
        "n_species": ("taxon_key", "nunique"),
        "n_detections_total": ("n_detections", "sum"),
        "B_H_index": ("biomass_contrib", "sum"),
        "M_H_index": ("metabolism_contrib", "sum"),
        "B_H_gt50": ("biomass_contrib_gt50", "sum"),
        "B_H_gt100": ("biomass_contrib_gt100", "sum"),
        "B_H_gt1000": ("biomass_contrib_gt1000", "sum"),
        "dominant_species": ("common_name", lambda x: str(x.dropna().mode().iloc[0]) if len(x.dropna().mode()) > 0 else ""),
        "species_list": ("common_name", lambda x: "; ".join(sorted(x.dropna().astype(str).unique()))),
    }

    site = (usable
            .groupby(["project_id", "deployment_id"])
            .agg(**agg_dict)
            .reset_index())

    site["megafauna_fraction"] = np.where(site["B_H_index"] > 0, (site["B_H_gt50"] / site["B_H_index"]) * 100, 0.0)

    meta_cols = ["project_id", "deployment_id", "longitude", "latitude", "trap_days"]
    if "project_name" in det.columns:
        meta_cols.append("project_name")
    
    all_deps = det[meta_cols].drop_duplicates()
    site = all_deps.merge(site, on=["project_id", "deployment_id"], how="left")

    site["n_species"] = site["n_species"].fillna(0).astype(int)
    return site


def process_region(region_dir: Path, region_name: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Iterate packages under a region and return clean detection and event frames."""
    print(f"\nProcessing region {region_name}...")
    packages = discover_wi_packages(region_dir)
    print(f"  Discovered {len(packages)} packages")
    
    all_deployments = []
    all_events = []
    
    for idx, p_dir in enumerate(packages):
        print(f"    [{idx+1}/{len(packages)}] Package: {p_dir.name}")
        try:
            deps = load_deployments(p_dir)
            evts = load_images(p_dir)
            all_deployments.append(deps)
            all_events.append(evts)
        except Exception as e:
            print(f"      ERROR loading package {p_dir.name}: {e}")
            
    if not all_deployments:
        return pd.DataFrame(), pd.DataFrame()
        
    deployments = pd.concat(all_deployments, ignore_index=True).drop_duplicates(subset=["project_id", "deployment_id"])
    events = pd.concat(all_events, ignore_index=True)
    dedup_cols = ["project_id", "deployment_id", "genus", "species", "timestamp", "number_of_objects"]
    dedup_cols = [c for c in dedup_cols if c in events.columns]
    events = events.drop_duplicates(subset=dedup_cols)

    print(f"  Total {region_name} clean deployments: {len(deployments)}")
    print(f"  Total {region_name} clean events:      {len(events)}")

    det = compute_detection_rates(deployments, events)
    det["region"] = region_name
    return det, events


# ── Spatial Clustering, Aggregation & GIS Preprocessing ────────────────────

def buffer_in_meters(geom, distance_meters):
    """Buffer a shapely geometry in meters using a local AEQD projection to avoid distortion."""
    centroid = geom.centroid
    lon, lat = centroid.x, centroid.y
    aeqd_proj = f"+proj=aeqd +lat_0={lat} +lon_0={lon} +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
    to_aeqd = Transformer.from_crs("EPSG:4326", aeqd_proj, always_xy=True)
    from_aeqd = Transformer.from_crs(aeqd_proj, "EPSG:4326", always_xy=True)
    geom_projected = transform(to_aeqd.transform, geom)
    geom_buffered_projected = geom_projected.buffer(distance_meters)
    return transform(from_aeqd.transform, geom_buffered_projected)


def haversine(lon1, lat1, lon2, lat2):
    """Calculate Haversine distance between two coordinates in km."""
    r = 6371.0
    rad = np.pi / 180
    dlon = (lon2 - lon1) * rad
    dlat = (lat2 - lat1) * rad
    lat1_r = lat1 * rad
    lat2_r = lat2 * rad
    a = np.sin(dlat / 2)**2 + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon / 2)**2
    c = 2 * np.arcsin(np.sqrt(a))
    return r * c


def get_master_cluster_map(df, threshold_km=None):
    """Create coordinate-to-cluster_id map using scipy single-linkage clustering."""
    if threshold_km is None:
        threshold_km = CLUSTER_THRESHOLD_KM
    unique_coords = df[['region', 'longitude', 'latitude']].drop_duplicates().reset_index(drop=True)
    unique_coords['cluster_id'] = ""
    
    for reg in unique_coords['region'].unique():
        sub = unique_coords[unique_coords['region'] == reg].copy()
        n = len(sub)
        if n == 0:
            continue
        
        if n > 1:
            coords = sub[['longitude', 'latitude']].values
            dist_matrix = np.zeros((n, n))
            for i in range(n):
                for j in range(i+1, n):
                    d = haversine(coords[i,0], coords[i,1], coords[j,0], coords[j,1])
                    dist_matrix[i,j] = d
                    dist_matrix[j,i] = d
            condensed = dist_matrix[np.triu_indices(n, k=1)]
            Z = linkage(condensed, method='single')
            labels = fcluster(Z, threshold_km, criterion='distance')
        else:
            labels = np.array([1])
            
        sub['c_num'] = labels
        sub['cluster_id'] = reg + "_" + sub['c_num'].astype(str).str.zfill(2)
        unique_coords.loc[unique_coords['region'] == reg, 'cluster_id'] = sub['cluster_id'].values
        
    return unique_coords.set_index(['region', 'longitude', 'latitude'])['cluster_id'].to_dict()


def aggregate_to_clusters(det, cluster_map):
    """Aggregate detections to cluster level and calculate biomass indices."""
    det = det.copy()
    det['cluster_id'] = det.apply(
        lambda r: cluster_map.get((r['region'], r['longitude'], r['latitude']), ""),
        axis=1
    )
    det = det[det['cluster_id'] != ""].copy()
    
    # Calculate total trap days per cluster
    deps = det[['region', 'cluster_id', 'project_id', 'deployment_id', 'trap_days']].drop_duplicates()
    cluster_trap_days = deps.groupby(['region', 'cluster_id'])['trap_days'].sum().to_dict()
    
    # Group detections by cluster and taxon key
    grp_cols = ["region", "cluster_id", "taxon_key", "class", "order", "family", "genus", "species", "taxon_quality"]
    agg_dict = {"n_detections": "sum"}
    for col in ["body_mass_kg", "common_name"]:
        if col in det.columns:
            agg_dict[col] = "first"
            
    clustered_det = det.groupby(grp_cols, dropna=False).agg(agg_dict).reset_index()
    clustered_det['cluster_trap_days'] = clustered_det.apply(
        lambda r: cluster_trap_days.get((r['region'], r['cluster_id']), 1.0),
        axis=1
    )
    
    # Cluster RAI
    clustered_det['RAI'] = (clustered_det['n_detections'] / clustered_det['cluster_trap_days']) * 100
    
    # Biophysical estimations
    usable = clustered_det[
        (clustered_det["taxon_quality"].isin(["species", "genus", "family"])) &
        (~clustered_det["taxon_quality"].isin(["blank", "human", "domestic"])) &
        (clustered_det["body_mass_kg"].notna()) &
        (clustered_det["body_mass_kg"] > 0)
    ].copy()
    
    usable["day_range_km"] = DAY_RANGE_COEFF * (usable["body_mass_kg"] ** DAY_RANGE_EXP)
    usable["corrected_RAI"] = usable["RAI"] / usable["day_range_km"]
    
    usable["biomass_contrib"] = usable["corrected_RAI"] * usable["body_mass_kg"]
    usable["metabolism_contrib"] = usable["corrected_RAI"] * FMR_COEFF * (usable["body_mass_kg"] ** FMR_EXP)
    
    usable["biomass_contrib_gt50"] = np.where(usable["body_mass_kg"] > 50.0, usable["biomass_contrib"], 0.0)
    usable["biomass_contrib_gt100"] = np.where(usable["body_mass_kg"] > 100.0, usable["biomass_contrib"], 0.0)
    usable["biomass_contrib_gt1000"] = np.where(usable["body_mass_kg"] > 1000.0, usable["biomass_contrib"], 0.0)
    
    agg_dict_metric = {
        "n_species": ("taxon_key", "nunique"),
        "n_detections_total": ("n_detections", "sum"),
        "B_H_index": ("biomass_contrib", "sum"),
        "M_H_index": ("metabolism_contrib", "sum"),
        "B_H_gt50": ("biomass_contrib_gt50", "sum"),
        "B_H_gt100": ("biomass_contrib_gt100", "sum"),
        "B_H_gt1000": ("biomass_contrib_gt1000", "sum"),
    }
    
    cluster_metrics = (usable
                       .groupby(["region", "cluster_id"])
                       .agg(**agg_dict_metric)
                       .reset_index())
    
    cluster_metrics["megafauna_fraction"] = np.where(
        cluster_metrics["B_H_index"] > 0,
        (cluster_metrics["B_H_gt50"] / cluster_metrics["B_H_index"]) * 100,
        0.0
    )
    # CR-5: Add new megafauna fraction metrics with appropriate thresholds
    cluster_metrics["megafauna_fraction_gt50"] = np.where(
        cluster_metrics["B_H_index"] > 0,
        (cluster_metrics["B_H_gt50"] / cluster_metrics["B_H_index"]) * 100,
        0.0
    )
    cluster_metrics["megafauna_fraction_gt100"] = np.where(
        cluster_metrics["B_H_index"] > 0,
        (cluster_metrics["B_H_gt100"] / cluster_metrics["B_H_index"]) * 100,
        0.0
    )
    
    # Calculate centroids and merge
    unique_dep_coords = det[['region', 'cluster_id', 'longitude', 'latitude']].drop_duplicates()
    centroids = unique_dep_coords.groupby(['region', 'cluster_id']).agg(
        longitude=('longitude', 'mean'),
        latitude=('latitude', 'mean')
    ).reset_index()
    
    all_meta = deps.groupby(['region', 'cluster_id']).agg(trap_days=('trap_days', 'sum')).reset_index()
    all_meta = all_meta.merge(centroids, on=['region', 'cluster_id'], how='left')
    
    cluster_metrics = all_meta.merge(cluster_metrics, on=["region", "cluster_id"], how="left")
    cluster_metrics["n_species"] = cluster_metrics["n_species"].fillna(0).astype(int)
    cluster_metrics["n_detections_total"] = cluster_metrics["n_detections_total"].fillna(0).astype(int)
    cluster_metrics["B_H_index"] = cluster_metrics["B_H_index"].fillna(0.0)
    cluster_metrics["M_H_index"] = cluster_metrics["M_H_index"].fillna(0.0)
    cluster_metrics["B_H_gt50"] = cluster_metrics["B_H_gt50"].fillna(0.0)
    cluster_metrics["B_H_gt100"] = cluster_metrics["B_H_gt100"].fillna(0.0)
    cluster_metrics["B_H_gt1000"] = cluster_metrics["B_H_gt1000"].fillna(0.0)
    
    return cluster_metrics, clustered_det


def check_gedi_5km_overlap(region, buffered_polygon):
    """Check if a buffered polygon overlaps at least one valid GEDI 5km pixel."""
    raster_path = OUTPUT_DIR / "EOdata" / f"analysis_stack_5000_{region}.tif"
    
    if not raster_path.exists():
        raise FileNotFoundError(
            f"Critical GEDI raster stack for {region} is missing! "
            f"Expected at outputs/EOdata/analysis_stack_5000_{region}.tif"
        )
        
    try:
        with rasterio.open(raster_path) as src:
            out_image, _ = mask(src, [buffered_polygon], crop=True, filled=False, indexes=3, all_touched=False)
            if isinstance(out_image, np.ma.MaskedArray):
                valid_pixels = np.sum(~out_image.mask & (out_image.data > 0))
            else:
                valid_pixels = np.sum(out_image > 0)
            return valid_pixels > 0
    except Exception as e:
        print(f"Warning masking GEDI raster: {e}")
        return False


def calculate_temporal_weights_py(det, cluster_map):
    """Calculate temporal weight (w_temp_cluster) for each cluster in Python."""
    det = det.copy()
    det['cluster_id'] = det.apply(
        lambda r: cluster_map.get((r['region'], r['longitude'], r['latitude']), ""),
        axis=1
    )
    det = det[det['cluster_id'] != ""].copy()
    
    deployments = det[['region', 'cluster_id', 'project_id', 'deployment_id', 'start_date', 'end_date', 'trap_days']].drop_duplicates().copy()
    deployments["start_date"] = pd.to_datetime(deployments["start_date"])
    gedi_start = pd.to_datetime("2019-04-17")
    
    deployments["years_before_gedi"] = (gedi_start - deployments["start_date"]).dt.days / 365.25
    
    brackets = CONFIG["temporal_decay"]["brackets"]
    weights = CONFIG["temporal_decay"]["weights"]
    
    def get_w_temp(yrs):
        if yrs <= brackets[0]:
            return weights[0]
        elif yrs <= brackets[1]:
            return weights[1]
        elif yrs <= brackets[2]:
            return weights[2]
        else:
            return weights[3]
            
    deployments["w_temp"] = deployments["years_before_gedi"].apply(get_w_temp)
    
    cluster_weights = deployments.groupby("cluster_id").apply(
        lambda g: (g["trap_days"] * g["w_temp"]).sum() / g["trap_days"].sum() if g["trap_days"].sum() > 0 else 1.0
    ).to_dict()
    
    return cluster_weights


def calculate_taxonomic_keep_proportions_py(det, cluster_map):
    """Calculate taxonomic keep proportion (p_keep) for each cluster in Python."""
    det = det.copy()
    det['cluster_id'] = det.apply(
        lambda r: cluster_map.get((r['region'], r['longitude'], r['latitude']), ""),
        axis=1
    )
    det = det[det['cluster_id'] != ""].copy()
    
    wild_det = det[~det["taxon_quality"].isin(["blank", "human", "domestic"])].copy()
    
    if len(wild_det) > 0:
        total_wild = wild_det.groupby("cluster_id")["n_detections"].sum()
        kept_wild = wild_det[wild_det["taxon_quality"].isin(["species", "genus", "family"])].groupby("cluster_id")["n_detections"].sum()
        p_keep_dict = (kept_wild / total_wild).fillna(1.0).to_dict()
    else:
        p_keep_dict = {}
        
    return p_keep_dict


# ── Main Script Entry ────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Process Congo, Amazon, and SE Asia camera trapping datasets with spatial clustering."
    )
    parser.add_argument("--congo-dir", type=Path, default=CONGO_DIR)
    parser.add_argument("--amazon-dir", type=Path, default=AMAZON_DIR)
    parser.add_argument("--seasia-dir", type=Path, default=SEASIA_DIR)
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

    # Process Southeast Asia
    seasia_det, seasia_raw_events = process_region(args.seasia_dir, "SE_Asia")

    # Combine detections
    print("\nMerging datasets...")
    joint_det = pd.concat([congo_det, amazon_det, seasia_det], ignore_index=True)

    # Match body masses
    print("Matching body masses...")
    joint_det = match_body_mass(joint_det, traits)
    
    match_summary = joint_det["mass_match_level"].value_counts()
    print("Match quality:")
    for lvl, cnt in match_summary.items():
        print(f"  {lvl:18s}: {cnt:5d} rows")

    # Run deployment-level site metrics (legacy metrics)
    print("\nComputing site-level heterotroph biomass and metabolism indices...")
    joint_metrics = compute_site_heterotroph_metrics(joint_det)
    joint_metrics["region"] = joint_metrics["project_id"].map(
        joint_det.drop_duplicates("project_id").set_index("project_id")["region"]
    )
    joint_metrics.to_csv(args.out_metrics, index=False)
    print(f"Saved legacy joint metrics (deployment-level) to {args.out_metrics}")

    # Spatial clustering and GEDI filtering
    if SPATIAL_LIBS_AVAILABLE:
        print(f"\nRunning single-linkage spatial clustering at {CLUSTER_THRESHOLD_KM} km on deployments...")
        cluster_map = get_master_cluster_map(joint_det, CLUSTER_THRESHOLD_KM)
        
        # Add cluster_id to event-level detections
        joint_det['cluster_id'] = joint_det.apply(
            lambda r: cluster_map.get((r['region'], r['longitude'], r['latitude']), ""),
            axis=1
        )
        # Save event-level detections with cluster assignments
        joint_det.to_csv(args.out_detections, index=False)
        print(f"Saved joint detections with cluster assignments to {args.out_detections}")

        print("Aggregating detections to mathematical cluster-level metrics...")
        cluster_metrics, clustered_det = aggregate_to_clusters(joint_det, cluster_map)

        print("\nCalculating temporal weights and keep proportions in Python...")
        p_keep_dict = calculate_taxonomic_keep_proportions_py(joint_det, cluster_map)
        w_temp_dict = calculate_temporal_weights_py(joint_det, cluster_map)

        cluster_metrics["p_keep"] = cluster_metrics["cluster_id"].map(p_keep_dict).fillna(1.0)
        cluster_metrics["w_temp_cluster"] = cluster_metrics["cluster_id"].map(w_temp_dict).fillna(1.0)

        print("\nFiltering clusters to only keep those overlapping at least one valid GEDI 5km pixel...")
        # Map deployments to cluster ids for grouping
        det_copy = joint_det.copy()
        valid_clusters = []
        for _, row in cluster_metrics.iterrows():
            c_id = row["cluster_id"]
            region = row["region"]
            
            c_deps = det_copy[det_copy["cluster_id"] == c_id]
            points = list(zip(c_deps["longitude"], c_deps["latitude"]))
            
            mp = MultiPoint(points)
            hull = mp.convex_hull
            buffered = buffer_in_meters(hull, BUFFER_METERS)
            
            if check_gedi_5km_overlap(region, buffered):
                valid_clusters.append(c_id)
                
        print(f"  ✓ Retained {len(valid_clusters)} / {len(cluster_metrics)} clusters with valid GEDI overlap.")
        cluster_metrics_filt = cluster_metrics[cluster_metrics["cluster_id"].isin(valid_clusters)].copy()
        
        # Save the aggregated cluster metrics
        cluster_metrics_path = OUTPUT_DIR / "camera_traps_cluster_level_metrics.csv"
        cluster_metrics_filt.to_csv(cluster_metrics_path, index=False)
        print(f"Saved aggregated cluster-level metrics to {cluster_metrics_path}")

        # Save robust metrics (trap_days >= MIN_TRAP_DAYS)
        robust_cluster_metrics = cluster_metrics_filt[cluster_metrics_filt["trap_days"] >= MIN_TRAP_DAYS].copy()
        robust_metrics_path = OUTPUT_DIR / "camera_traps_cluster_level_metrics_robust.csv"
        robust_cluster_metrics.to_csv(robust_metrics_path, index=False)
        print(f"Saved robust cluster-level metrics to {robust_metrics_path}")

        # Generate robust cluster buffered MCPs as GeoJSON
        print("Generating and saving robust cluster buffered MCPs as GeoJSON...")
        robust_det_for_geojson = det_copy[det_copy['cluster_id'].isin(robust_cluster_metrics['cluster_id'])].copy()
        robust_det_for_geojson.to_csv(OUTPUT_DIR / "camera_traps_robust_detections.csv", index=False)
        print(f"Saved robust filtered detections to {OUTPUT_DIR / 'camera_traps_robust_detections.csv'}")

        # Group to unique deployments for polygons
        deps_unique = (robust_det_for_geojson.groupby(["region", "cluster_id", "project_id", "deployment_id"])
                       .agg(lon=("longitude", "first"),
                            lat=("latitude", "first"),
                            trap_days=("trap_days", "first"))
                       .reset_index())

        features = []
        for c_id in robust_cluster_metrics['cluster_id'].unique():
            c_deps = deps_unique[deps_unique["cluster_id"] == c_id]
            if len(c_deps) == 0:
                continue
            points = list(zip(c_deps["lon"], c_deps["lat"]))

            mp = MultiPoint(points)
            hull = mp.convex_hull
            buffered = buffer_in_meters(hull, BUFFER_METERS)

            c_info = robust_cluster_metrics[robust_cluster_metrics['cluster_id'] == c_id].iloc[0]

            properties = {
                "cluster_id": str(c_id),
                "region": str(c_info["region"]),
                "trap_days": float(c_info["trap_days"]),
                "n_species": int(c_info["n_species"]),
                "B_H_index": float(c_info["B_H_index"]),
                "M_H_index": float(c_info["M_H_index"]),
                "B_H_gt50": float(c_info["B_H_gt50"]),
                "B_H_gt100": float(c_info["B_H_gt100"]),
                "B_H_gt1000": float(c_info["B_H_gt1000"]),
                "megafauna_fraction": float(c_info["megafauna_fraction"]),
                "megafauna_fraction_gt50": float(c_info["megafauna_fraction_gt50"]),
                "megafauna_fraction_gt100": float(c_info["megafauna_fraction_gt100"]),
                "p_keep": float(c_info["p_keep"]),
                "w_temp_cluster": float(c_info["w_temp_cluster"])
            }

            feature = {
                "type": "Feature",
                "geometry": mapping(buffered),
                "properties": properties
            }
            features.append(feature)

        geojson = {
            "type": "FeatureCollection",
            "features": features
        }

        geojson_path = OUTPUT_DIR / "camera_traps_robust_buffered_mcps.geojson"
        with open(geojson_path, "w") as f:
            json.dump(geojson, f, indent=2)
        print(f"Saved robust cluster buffered MCPs to {geojson_path}")

    else:
        # Fallback if spatial dependencies are missing
        print("\nWARNING: scipy, shapely, or rasterio are not available.")
        print("  Saving raw detections only, spatial clustering was skipped.")
        joint_det.to_csv(args.out_detections, index=False)


if __name__ == "__main__":
    main()
