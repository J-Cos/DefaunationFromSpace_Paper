#!/usr/bin/env python3
"""
visualise_camera_traps.py
=========================
Generates PNAS-style publication-quality multi-panel figures comparing
vertebrate community structure and biophysical scaling across the
Amazon and Congo Basins. All figures use 11.1 km single-linkage spatial
clusters as independent units rather than individual deployments.

Figure 1: Camera Trap Detections & Diversity Comparison
    (A) Maps of fuzzed camera trap MCPs (+ 11.1km buffer) colored by cumulative trap days
    (B) Rank-abundance curves of species in both basins
    (C) Taxonomic composition by order (proportional horizontal bars)
    (D) Cluster-level wild heterotroph RAI distributions (box/violin plot)

Figure 2: Biophysical Scaling & Energetics Comparison (Cluster-level)
    (A) log B_H vs log M_H scatter plot with fitted scaling lines
    (B) Species richness vs log B_H scatter plot
    (C) Standing heterotroph biomass (B_H) index distribution
    (D) Standing heterotroph metabolism (M_H) index distribution

Figure 3: Vertebrate Body Size & Megafaunal Biomass Comparison (Cluster-level)
    (A) standing Biomass Index of wild animals > 50 kg per cluster
    (B) standing Biomass Index of megafauna > 100 kg per cluster
    (C) Log body mass probability density curves of detected independent events
    (D) Megafauna biomass fraction (%) of total standing biomass

Output: - figures/camera_traps_fig01_detections.pdf/png
        - figures/camera_traps_fig02_biophysical.pdf/png
        - figures/camera_traps_fig03_bodysize.pdf/png
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.colors as colors
from pathlib import Path
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, fcluster
from shapely.geometry import MultiPoint, shape
import rasterio
from rasterio.mask import mask
import json

# ── Paths ───────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
OUTPUT_DIR = PROJECT_DIR / "outputs"
FIG_DIR = PROJECT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ── PNAS style configuration ───────────────────────────────────────────────
SINGLE_COL = 3.42
DOUBLE_COL = 7.0

PAL = {
    "Congo":   "#2166AC",  # Blue for Congo
    "Amazon":  "#B2182B",  # Red for Amazon
    "SE_Asia": "#1B7837",  # Green for SE Asia
    "green":   "#1B7837",
    "orange":  "#E08214",
    "purple":  "#6A3D9A",
    "grey":    "#878787",
    "teal":    "#35978F",
    "pink":    "#C51B7D",
    "light_grey": "#E0E0E0",
}

ORDER_COLOURS = {
    "Cetartiodactyla":  "#1B7837",  # Green
    "Proboscidea":      "#2166AC",  # Blue
    "Carnivora":        "#E08214",  # Orange
    "Rodentia":         "#6A3D9A",  # Purple
    "Primates":         "#C51B7D",  # Pink
    "Cingulata":        "#35978F",  # Teal
    "Pilosa":           "#DFC27D",  # Gold
    "Didelphimorphia":  "#8c510a",  # Brown
    "Galliformes":      "#543005",
    "Columbiformes":    "#b2182b",
    "Perissodactyla":   "#4d4d4d",
    "Other":            "#878787",  # Grey
}


def set_pnas_style():
    """Apply PNAS-compatible matplotlib rcParams."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 7,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,
        "axes.linewidth": 0.5,
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size": 3,
        "ytick.major.size": 3,
        "xtick.minor.size": 1.5,
        "ytick.minor.size": 1.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "patch.linewidth": 0.5,
    })


def panel_label(ax, label, x=-0.12, y=1.05):
    """Add bold panel label (A, B, C...) to an axes."""
    ax.text(x, y, label, transform=ax.transAxes,
            fontsize=9, fontweight="bold", va="top", ha="left")


# ── Coordinate mapping helper ────────────────────────────────────────────────

def get_coordinate_cluster_map(df):
    """Build a mapping from (region, lon, lat) to cluster_id directly from data."""
    return df.drop_duplicates(['region', 'longitude', 'latitude']).set_index(['region', 'longitude', 'latitude'])['cluster_id'].to_dict()





def _compute_region_extent(det: pd.DataFrame, region: str, pad: float = 1.0) -> list:
    """Dynamically compute [xmin, xmax, ymin, ymax] map extent from deployment coords."""
    sub = det[det["region"] == region]
    if sub.empty:
        return [0, 1, 0, 1]
    return [
        sub["longitude"].min() - pad,
        sub["longitude"].max() + pad,
        sub["latitude"].min() - pad,
        sub["latitude"].max() + pad
    ]


# ── Figure 1: Camera Trap Detections & Diversity ───────────────────────────

def make_figure1(det: pd.DataFrame, cluster_metrics: pd.DataFrame, geojson_features: list, fig_dir: Path):
    """
    4-panel summary of camera trap detections and diversity across basins.
    (A) Deployment maps using Minimum Convex Polygons + 11.1km buffer colored by effort
    (B) Rank-abundance comparison
    (C) Proportional taxonomic order composition
    (D) Cluster-level RAI distribution (landscape scale)
    """
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        cartopy_available = True
    except ImportError:
        cartopy_available = False

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.98))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.05], hspace=0.45, wspace=0.35)

    # ── (A) Geographic Maps (MCPs + 11.1km Buffers colored by trap_days) ────
    gs_maps = gs[0, :].subgridspec(1, 3, wspace=0.15)
    
    regions_info = [
        {"name": "Congo", "extent": _compute_region_extent(det, "Congo"), "gs": gs_maps[0]},
        {"name": "Amazon", "extent": _compute_region_extent(det, "Amazon"), "gs": gs_maps[1]},
        {"name": "SE_Asia", "extent": _compute_region_extent(det, "SE_Asia"), "gs": gs_maps[2]}
    ]

    # Shared log-scaled effort normalization across all clusters globally
    min_days = max(10, cluster_metrics['trap_days'].min())
    max_days = cluster_metrics['trap_days'].max()
    norm = colors.LogNorm(vmin=min_days, vmax=max_days)
    cmap = plt.cm.viridis

    for idx, reg in enumerate(regions_info):
        if cartopy_available:
            ax_map = fig.add_subplot(reg["gs"], projection=ccrs.PlateCarree())
            ax_map.set_extent(reg["extent"], crs=ccrs.PlateCarree())
            
            ax_map.add_feature(cfeature.LAND, facecolor="#F8F8F6", zorder=0)
            ax_map.add_feature(cfeature.OCEAN, facecolor="#EBF4F6", zorder=0)
            ax_map.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor="#999999", zorder=1)
            ax_map.add_feature(cfeature.BORDERS, linewidth=0.2, linestyle=":", edgecolor="#B0B0B0", zorder=1)
            ax_map.add_feature(cfeature.RIVERS, linewidth=0.2, edgecolor="#D0E0E3", zorder=1)
            
            gl = ax_map.gridlines(draw_labels=True, linewidth=0.1, color="#DDDDDD", alpha=0.5)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {"size": 4.5}
            gl.ylabel_style = {"size": 4.5}
        else:
            ax_map = fig.add_subplot(reg["gs"])
            ax_map.set_xlim(reg["extent"][0], reg["extent"][1])
            ax_map.set_ylim(reg["extent"][2], reg["extent"][3])
            ax_map.set_facecolor("#FAFAFA")
            ax_map.set_xlabel("Longitude (°E)", fontsize=5.5)
            ax_map.set_ylabel("Latitude (°N)", fontsize=5.5)
            ax_map.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0")

        # Plot pre-computed projection-safe buffered MCPs from GeoJSON
        n_clusters = 0
        for feat in geojson_features:
            if feat["region"] != reg["name"]:
                continue
            n_clusters += 1
            
            buffered = feat["geometry"]
            c_days = feat["trap_days"]
            color = cmap(norm(c_days))
            
            try:
                if cartopy_available:
                    ax_map.add_geometries(
                        [buffered], crs=ccrs.PlateCarree(),
                        facecolor=color, edgecolor="black", linewidth=0.3,
                        alpha=0.6, zorder=3
                    )
                else:
                    if buffered.geom_type == 'Polygon':
                        x, y = buffered.exterior.xy
                        ax_map.fill(x, y, facecolor=color, edgecolor="black", linewidth=0.3, alpha=0.6, zorder=3)
                    elif buffered.geom_type == 'MultiPolygon':
                        for poly in buffered.geoms:
                            x, y = poly.exterior.xy
                            ax_map.fill(x, y, facecolor=color, edgecolor="black", linewidth=0.3, alpha=0.6, zorder=3)
            except Exception as e:
                print(f"      Warning: Drawing geometry failed: {e}")

        ax_map.set_title(f"{reg['name']} Basin ({n_clusters} clusters)", fontsize=7.5, pad=3)
        if idx == 0:
            panel_label(ax_map, "A", x=-0.08, y=1.05)

    # shared horizontal colorbar for map effort below top row
    cbar_ax = fig.add_axes([0.35, 0.54, 0.30, 0.012])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                      cax=cbar_ax, orientation='horizontal')
    cb.set_label('Cumulative effort (trap-days)', fontsize=5.5)
    cb.ax.tick_params(labelsize=4.5)

    # ── (B) Rank-Abundance Curves ───────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub_det = det[(det["region"] == region) & (det["taxon_quality"] == "species")].copy()
        
        sp_counts = (sub_det.groupby(["genus", "species", "common_name"])["n_detections"]
                     .sum().sort_values(ascending=False).reset_index())
        sp_counts["rank"] = range(1, len(sp_counts) + 1)
        
        ax.plot(sp_counts["rank"], sp_counts["n_detections"], 
                color=PAL[region], lw=1.2, marker="o", ms=2.5, alpha=0.8,
                label=f"{region} (S={len(sp_counts)})")
        
        top_spp = sp_counts.iloc[0]
        ax.annotate(top_spp["common_name"], 
                    xy=(1, top_spp["n_detections"]), 
                    xytext=(8, -2), textcoords="offset points",
                    fontsize=4.5, fontstyle="italic", ha="left", va="center",
                    arrowprops=dict(arrowstyle="->", lw=0.3, color=PAL[region]))
        
    ax.set_xlabel("Species rank")
    ax.set_ylabel("Total detections")
    ax.set_yscale("log")
    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    ax.legend(frameon=False, loc="upper right")
    panel_label(ax, "B")

    # ── (C) Taxonomic Order Composition ──────────────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    order_data = []
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub_det = det[(det["region"] == region) & (det["taxon_quality"].isin(["species", "genus", "family"]))].copy()
        sub_det = sub_det[~sub_det["taxon_quality"].isin(["blank", "human", "domestic"])]
        
        counts = sub_det.groupby("order")["n_detections"].sum()
        total = counts.sum()
        pct = (counts / total * 100).to_dict()
        
        for ord_name, ord_pct in pct.items():
            order_data.append({
                "region": region,
                "order": ord_name if ord_name in ORDER_COLOURS else "Other",
                "pct": ord_pct
            })
            
    df_ord = pd.DataFrame(order_data)
    df_pivot = df_ord.pivot_table(index="region", columns="order", values="pct", aggfunc="sum").fillna(0)
    
    legend_orders = [o for o in ORDER_COLOURS if o in df_pivot.columns]
    if "Other" in df_pivot.columns and "Other" not in legend_orders:
        legend_orders.append("Other")
    df_pivot = df_pivot[legend_orders]

    left = np.zeros(3)
    regions_idx = ["Congo", "Amazon", "SE_Asia"]
    
    for ord_name in legend_orders:
        pcts = df_pivot[ord_name].reindex(regions_idx).values
        ax.barh(regions_idx, pcts, left=left, color=ORDER_COLOURS[ord_name],
                edgecolor="white", linewidth=0.3, height=0.55, label=ord_name)
        left += pcts
        
    ax.set_xlabel("Proportion of detections (%)")
    ax.set_xlim(0, 100)
    ax.set_title("Vertebrate order composition", fontsize=7.5, pad=3)
    ax.legend(bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=False, 
              handlelength=0.8, handleheight=0.8, borderpad=0.1)
    panel_label(ax, "C")

    # ── (D) Cluster-level Wild Heterotroph RAI Distribution ──────────────────
    ax = fig.add_subplot(gs[1, 2])
    
    # Calculate wild RAI at the cluster level
    wild_clustered = cluster_metrics[['region', 'cluster_id', 'n_detections_total', 'trap_days']].copy()
    wild_clustered['RAI'] = (wild_clustered['n_detections_total'] / wild_clustered['trap_days']) * 100
    
    box_data = [
        wild_clustered[wild_clustered["region"] == "Congo"]["RAI"].values,
        wild_clustered[wild_clustered["region"] == "Amazon"]["RAI"].values,
        wild_clustered[wild_clustered["region"] == "SE_Asia"]["RAI"].values
    ]
    
    bp = ax.boxplot(box_data, tick_labels=["Congo", "Amazon", "SE_Asia"], patch_artist=True,
                    widths=0.45, showfliers=False, zorder=2)
    
    colors_bp = [PAL["Congo"], PAL["Amazon"], PAL["SE_Asia"]]
    for patch, color in zip(bp["boxes"], colors_bp):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)
        
    for element in ["whiskers", "caps", "medians"]:
        plt.setp(bp[element], color="black", lw=0.6)
        
    # Jitter raw clusters
    for i, region in enumerate(["Congo", "Amazon", "SE_Asia"]):
        vals = box_data[i]
        x_jitter = np.random.default_rng(i).normal(i + 1, 0.04, len(vals))
        ax.scatter(x_jitter, vals, color=colors_bp[i], s=5, alpha=0.6, edgecolors="none", zorder=3)
        
        median_val = np.median(vals)
        ax.text(i + 1, median_val * 1.15 if median_val > 0 else 5, f"med={median_val:.1f}",
                ha="center", va="bottom", fontsize=5, color="black", fontweight="bold")
        
    ax.set_ylabel("Cluster-level wild RAI\n(detections / 100 trap-days)")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5, zorder=1)
    panel_label(ax, "D")

    # Save
    for ext in ["pdf", "png"]:
        fig.savefig(fig_dir / f"camera_traps_fig01_detections.{ext}", dpi=300)
    plt.close(fig)
    print(f"Saved landscape cluster-level camera_traps_fig01_detections.pdf/png.")


# ── Figure 2: Biophysical Scaling & Energetics ─────────────────────────────

def make_figure2(valid: pd.DataFrame, fig_dir: Path):
    """
    4-panel biophysical scaling and energetics summary using mathematical clusters.
    (A) log B_H vs log M_H scatter plot (aggregated clusters)
    (B) Cluster-level species richness vs log B_H scatter plot
    (C) Standing heterotroph biomass (B_H) index distribution
    (D) Standing heterotroph metabolism (M_H) index distribution
    """
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.85), constrained_layout=True)

    # Filter out valid active clusters with non-zero biomass
    valid = valid[(valid["B_H_index"] > 0) & (valid["M_H_index"] > 0)].copy()

    # ── (A) log B_H vs log M_H Scatter ──────────────────────────────────────
    ax = axes[0, 0]
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub = valid[valid["region"] == region]
        ax.scatter(np.log10(sub["B_H_index"]), np.log10(sub["M_H_index"]),
                   c=PAL[region], s=16, alpha=0.7, edgecolors="white", linewidths=0.2,
                   label=region, zorder=3)
        
    # Fit pooled line
    logB = np.log10(valid["B_H_index"])
    logM = np.log10(valid["M_H_index"])
    slope, intercept = np.polyfit(logB, logM, 1)
    
    xfit = np.linspace(logB.min(), logB.max(), 50)
    yfit = slope * xfit + intercept
    ax.plot(xfit, yfit, color=PAL["grey"], ls="--", lw=0.9, zorder=2)
    
    corr = valid["B_H_index"].corr(valid["M_H_index"], method="spearman")
    ax.text(0.05, 0.92, f"slope = {slope:.3f}\n$r_S$ = {corr:.3f}", 
            transform=ax.transAxes, fontsize=5.5, color="black", fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1))
    
    ax.set_xlabel("log₁₀ standing Biomass Index ($B_H$)")
    ax.set_ylabel("log₁₀ standing Metabolism Index ($M_H$)")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    ax.legend(frameon=False, loc="lower right")
    panel_label(ax, "A")

    # ── (B) Cluster-level Species Richness vs log B_H ────────────────────────
    ax = axes[0, 1]
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub = valid[valid["region"] == region]
        ax.scatter(sub["n_species"], np.log10(sub["B_H_index"]),
                   c=PAL[region], s=16, alpha=0.7, edgecolors="white", linewidths=0.2,
                   zorder=3)
        
        if len(sub) > 2:
            s_slope, s_int = np.polyfit(sub["n_species"], np.log10(sub["B_H_index"]), 1)
            x_s = np.linspace(sub["n_species"].min(), sub["n_species"].max(), 50)
            ax.plot(x_s, s_slope * x_s + s_int, color=PAL[region], ls=":", lw=0.8, zorder=2)

    ax.set_xlabel("Cluster-level species richness ($S$)")
    ax.set_ylabel("log₁₀ standing Biomass Index ($B_H$)")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    panel_label(ax, "B")

    # ── (C) Standing Biomass (B_H) Distribution ────────────────────────────
    ax = axes[1, 0]
    logB_congo = np.log10(valid[valid["region"] == "Congo"]["B_H_index"])
    logB_amazon = np.log10(valid[valid["region"] == "Amazon"]["B_H_index"])
    logB_seasia = np.log10(valid[valid["region"] == "SE_Asia"]["B_H_index"])
    
    bins = np.histogram(np.hstack((logB_congo, logB_amazon, logB_seasia)), bins=12)[1]
    
    ax.hist(logB_congo, bins=bins, color=PAL["Congo"], alpha=0.45, edgecolor=PAL["Congo"],
            linewidth=0.4, label="Congo", density=False)
    ax.hist(logB_amazon, bins=bins, color=PAL["Amazon"], alpha=0.45, edgecolor=PAL["Amazon"],
            linewidth=0.4, label="Amazon", density=False)
    ax.hist(logB_seasia, bins=bins, color=PAL["SE_Asia"], alpha=0.45, edgecolor=PAL["SE_Asia"],
            linewidth=0.4, label="SE_Asia", density=False)
    
    med_c = np.median(logB_congo)
    med_a = np.median(logB_amazon)
    med_s = np.median(logB_seasia)
    
    ax.axvline(med_c, color=PAL["Congo"], ls="--", lw=0.8, zorder=4)
    ax.axvline(med_a, color=PAL["Amazon"], ls="--", lw=0.8, zorder=4)
    ax.axvline(med_s, color=PAL["SE_Asia"], ls="--", lw=0.8, zorder=4)
    
    ax.text(med_c + 0.05, ax.get_ylim()[1] * 0.85, f"med={10**med_c:.1f}", 
            color=PAL["Congo"], fontsize=5, fontweight="bold")
    ax.text(med_a - 0.05, ax.get_ylim()[1] * 0.70, f"med={10**med_a:.1f}", 
            color=PAL["Amazon"], fontsize=5, fontweight="bold", ha="right")
    ax.text(med_s + 0.05, ax.get_ylim()[1] * 0.55, f"med={10**med_s:.1f}", 
            color=PAL["SE_Asia"], fontsize=5, fontweight="bold")
    
    ax.set_xlabel("log₁₀ standing Biomass Index ($B_H$)")
    ax.set_ylabel("Number of spatial clusters")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    panel_label(ax, "C")

    # ── (D) Standing Metabolism (M_H) Distribution ─────────────────────────
    ax = axes[1, 1]
    logM_congo = np.log10(valid[valid["region"] == "Congo"]["M_H_index"])
    logM_amazon = np.log10(valid[valid["region"] == "Amazon"]["M_H_index"])
    logM_seasia = np.log10(valid[valid["region"] == "SE_Asia"]["M_H_index"])
    
    bins_m = np.histogram(np.hstack((logM_congo, logM_amazon, logM_seasia)), bins=12)[1]
    
    ax.hist(logM_congo, bins=bins_m, color=PAL["Congo"], alpha=0.45, edgecolor=PAL["Congo"],
            linewidth=0.4, label="Congo", density=False)
    ax.hist(logM_amazon, bins=bins_m, color=PAL["Amazon"], alpha=0.45, edgecolor=PAL["Amazon"],
            linewidth=0.4, label="Amazon", density=False)
    ax.hist(logM_seasia, bins=bins_m, color=PAL["SE_Asia"], alpha=0.45, edgecolor=PAL["SE_Asia"],
            linewidth=0.4, label="SE_Asia", density=False)
    
    med_mc = np.median(logM_congo)
    med_ma = np.median(logM_amazon)
    med_ms = np.median(logM_seasia)
    ax.axvline(med_mc, color=PAL["Congo"], ls="--", lw=0.8, zorder=4)
    ax.axvline(med_ma, color=PAL["Amazon"], ls="--", lw=0.8, zorder=4)
    ax.axvline(med_ms, color=PAL["SE_Asia"], ls="--", lw=0.8, zorder=4)
    
    ax.text(med_mc + 0.05, ax.get_ylim()[1] * 0.85, f"med={10**med_mc:.1f}", 
            color=PAL["Congo"], fontsize=5, fontweight="bold")
    ax.text(med_ma - 0.05, ax.get_ylim()[1] * 0.70, f"med={10**med_ma:.1f}", 
            color=PAL["Amazon"], fontsize=5, fontweight="bold", ha="right")
    ax.text(med_ms + 0.05, ax.get_ylim()[1] * 0.55, f"med={10**med_ms:.1f}", 
            color=PAL["SE_Asia"], fontsize=5, fontweight="bold")
    
    ax.set_xlabel("log₁₀ standing Metabolism Index ($M_H$)")
    ax.set_ylabel("Number of spatial clusters")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    panel_label(ax, "D")

    # Save
    for ext in ["pdf", "png"]:
        fig.savefig(fig_dir / f"camera_traps_fig02_biophysical.{ext}", dpi=300)
    plt.close(fig)
    print(f"Saved cluster-level camera_traps_fig02_biophysical.pdf/png.")


# ── Figure 3: Vertebrate Body Size & Megafauna ─────────────────────────────

def make_figure3(det: pd.DataFrame, cluster_metrics: pd.DataFrame, geojson_features: list, fig_dir: Path):
    """
    5-panel body size and megafaunal comparison (Landscape Cluster-level).
    (A) Geographic maps of megafauna biomass ($B_{H, >50}$) for Congo & Amazon clusters using buffered MCPs (0.05 degree)
    (B) standing Biomass Index of wild animals > 50 kg per cluster
    (C) standing Biomass Index of megafauna > 100 kg per cluster
    (D) Log body mass probability density curves of detected independent events
    (E) Megafauna biomass fraction (%) of total standing biomass
    """
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        cartopy_available = True
    except ImportError:
        cartopy_available = False

    # 5-panel layout
    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 1.35))
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.0], hspace=0.4, wspace=0.3)

    # ── (A) Geographic Maps (Buffered MCPs from GeoJSON colored by B_H_gt50) ────
    gs_maps = gs[0, :].subgridspec(1, 3, wspace=0.15)
    
    regions_info = [
        {"name": "Congo", "extent": _compute_region_extent(det, "Congo"), "gs": gs_maps[0]},
        {"name": "Amazon", "extent": _compute_region_extent(det, "Amazon"), "gs": gs_maps[1]},
        {"name": "SE_Asia", "extent": _compute_region_extent(det, "SE_Asia"), "gs": gs_maps[2]}
    ]

    # Shared log-scaled biomass normalization across all clusters globally
    biomass_map = {(r["region"], r["cluster_id"]): r["B_H_gt50"] for _, r in cluster_metrics.iterrows()}
    
    non_zero_biomass = cluster_metrics[cluster_metrics["B_H_gt50"] > 0]["B_H_gt50"]
    min_biomass = non_zero_biomass.min() if len(non_zero_biomass) > 0 else 1.0
    max_biomass = non_zero_biomass.max() if len(non_zero_biomass) > 0 else 100.0
    if min_biomass <= 0:
        min_biomass = 1.0
    
    norm = colors.LogNorm(vmin=min_biomass, vmax=max_biomass)
    cmap = plt.cm.plasma

    for idx, reg in enumerate(regions_info):
        if cartopy_available:
            ax_map = fig.add_subplot(reg["gs"], projection=ccrs.PlateCarree())
            ax_map.set_extent(reg["extent"], crs=ccrs.PlateCarree())
            
            ax_map.add_feature(cfeature.LAND, facecolor="#F8F8F6", zorder=0)
            ax_map.add_feature(cfeature.OCEAN, facecolor="#EBF4F6", zorder=0)
            ax_map.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor="#999999", zorder=1)
            ax_map.add_feature(cfeature.BORDERS, linewidth=0.2, linestyle=":", edgecolor="#B0B0B0", zorder=1)
            ax_map.add_feature(cfeature.RIVERS, linewidth=0.2, edgecolor="#D0E0E3", zorder=1)
            
            gl = ax_map.gridlines(draw_labels=True, linewidth=0.1, color="#DDDDDD", alpha=0.5)
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {"size": 4.5}
            gl.ylabel_style = {"size": 4.5}
        else:
            ax_map = fig.add_subplot(reg["gs"])
            ax_map.set_xlim(reg["extent"][0], reg["extent"][1])
            ax_map.set_ylim(reg["extent"][2], reg["extent"][3])
            ax_map.set_facecolor("#FAFAFA")
            ax_map.set_xlabel("Longitude (°E)", fontsize=5.5)
            ax_map.set_ylabel("Latitude (°N)", fontsize=5.5)
            ax_map.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0")

        # Plot pre-computed projection-safe buffered MCPs from GeoJSON
        n_clusters = 0
        for feat in geojson_features:
            if feat["region"] != reg["name"]:
                continue
            n_clusters += 1
            
            c_id = feat["cluster_id"]
            buffered = feat["geometry"]
            
            c_id_int = int(c_id) if isinstance(c_id, str) and c_id.isdigit() else c_id
            c_biomass = biomass_map.get((reg["name"], c_id), biomass_map.get((reg["name"], c_id_int), 0.0))
            
            if c_biomass > 0.0:
                color = cmap(norm(c_biomass))
            else:
                color = "#E0E0E0"  # light grey for absence
            
            try:
                if cartopy_available:
                    ax_map.add_geometries(
                        [buffered], crs=ccrs.PlateCarree(),
                        facecolor=color, edgecolor="black", linewidth=0.3,
                        alpha=0.6, zorder=3
                    )
                else:
                    if buffered.geom_type == 'Polygon':
                        x, y = buffered.exterior.xy
                        ax_map.fill(x, y, facecolor=color, edgecolor="black", linewidth=0.3, alpha=0.6, zorder=3)
                    elif buffered.geom_type == 'MultiPolygon':
                        for poly in buffered.geoms:
                            x, y = poly.exterior.xy
                            ax_map.fill(x, y, facecolor=color, edgecolor="black", linewidth=0.3, alpha=0.6, zorder=3)
            except Exception as e:
                print(f"      Warning: Shapely buffer failed for {c_id}: {e}")

        ax_map.set_title(f"{reg['name']} Basin ({len(unique_c_ids)} clusters)", fontsize=7.5, pad=3)
        if idx == 0:
            panel_label(ax_map, "A", x=-0.08, y=1.05)

    # shared horizontal colorbar for map biomass below top row
    cbar_ax = fig.add_axes([0.35, 0.655, 0.30, 0.012])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                      cax=cbar_ax, orientation='horizontal')
    cb.set_label('Megafauna biomass ($B_{H, >50}$, relative units)', fontsize=5.5)
    cb.ax.tick_params(labelsize=4.5)

    # Boxplot colors
    congo = cluster_metrics[cluster_metrics["region"] == "Congo"]
    amazon = cluster_metrics[cluster_metrics["region"] == "Amazon"]
    seasia = cluster_metrics[cluster_metrics["region"] == "SE_Asia"]
    colors_6 = [PAL["Congo"], "#B2D2EC", PAL["Amazon"], "#F4A582", PAL["SE_Asia"], "#A1D99B"]

    # ── (B) Biomass Index > 50 kg ───────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    congo_pres = congo[congo["B_H_gt50"] > 0]["B_H_gt50"].values
    congo_abs = congo[congo["B_H_gt50"] == 0]["B_H_gt50"].values
    amazon_pres = amazon[amazon["B_H_gt50"] > 0]["B_H_gt50"].values
    amazon_abs = amazon[amazon["B_H_gt50"] == 0]["B_H_gt50"].values
    seasia_pres = seasia[seasia["B_H_gt50"] > 0]["B_H_gt50"].values
    seasia_abs = seasia[seasia["B_H_gt50"] == 0]["B_H_gt50"].values
    
    box_data_50 = [congo_pres, congo_abs, amazon_pres, amazon_abs, seasia_pres, seasia_abs]
    labels_50 = [
        f"Congo\nPres\n(n={len(congo_pres)})",
        f"Congo\nAbs\n(n={len(congo_abs)})",
        f"Amazon\nPres\n(n={len(amazon_pres)})",
        f"Amazon\nAbs\n(n={len(amazon_abs)})",
        f"SE_Asia\nPres\n(n={len(seasia_pres)})",
        f"SE_Asia\nAbs\n(n={len(seasia_abs)})"
    ]
    
    bp_50 = ax.boxplot(box_data_50, tick_labels=labels_50, patch_artist=True,
                       widths=0.5, showfliers=False, zorder=2)
    
    for patch, color in zip(bp_50["boxes"], colors_6):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.5)
    for element in ["whiskers", "caps", "medians"]:
        plt.setp(bp_50[element], color="black", lw=0.5)
        
    for i, vals in enumerate(box_data_50):
        if len(vals) > 0:
            x_jitter = np.random.default_rng(i).normal(i + 1, 0.04, len(vals))
            ax.scatter(x_jitter, vals, color=colors_6[i], s=3.5, alpha=0.6, edgecolors="none", zorder=3)
            med_val = np.median(vals)
            ax.text(i + 1, med_val + 5 if med_val > 0 else 5, f"med={med_val:.1f}",
                    ha="center", va="bottom", fontsize=5.0, color="black", fontweight="bold")
        
    ax.set_ylabel("Biomass Index for animals >50 kg ($B_{H,>50}$, relative units)")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5, zorder=1)
    ax.tick_params(axis='x', labelsize=4.8)
    panel_label(ax, "B")

    # ── (C) Biomass Index > 100 kg ──────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 1])
    congo_pres_100 = congo[congo["B_H_gt100"] > 0]["B_H_gt100"].values
    congo_abs_100 = congo[congo["B_H_gt100"] == 0]["B_H_gt100"].values
    amazon_pres_100 = amazon[amazon["B_H_gt100"] > 0]["B_H_gt100"].values
    amazon_abs_100 = amazon[amazon["B_H_gt100"] == 0]["B_H_gt100"].values
    seasia_pres_100 = seasia[seasia["B_H_gt100"] > 0]["B_H_gt100"].values
    seasia_abs_100 = seasia[seasia["B_H_gt100"] == 0]["B_H_gt100"].values
    
    box_data_100 = [congo_pres_100, congo_abs_100, amazon_pres_100, amazon_abs_100, seasia_pres_100, seasia_abs_100]
    labels_100 = [
        f"Congo\nPres\n(n={len(congo_pres_100)})",
        f"Congo\nAbs\n(n={len(congo_abs_100)})",
        f"Amazon\nPres\n(n={len(amazon_pres_100)})",
        f"Amazon\nAbs\n(n={len(amazon_abs_100)})",
        f"SE_Asia\nPres\n(n={len(seasia_pres_100)})",
        f"SE_Asia\nAbs\n(n={len(seasia_abs_100)})"
    ]
    
    bp_100 = ax.boxplot(box_data_100, tick_labels=labels_100, patch_artist=True,
                        widths=0.5, showfliers=False, zorder=2)
    
    for patch, color in zip(bp_100["boxes"], colors_6):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.5)
    for element in ["whiskers", "caps", "medians"]:
        plt.setp(bp_100[element], color="black", lw=0.5)
        
    for i, vals in enumerate(box_data_100):
        if len(vals) > 0:
            x_jitter = np.random.default_rng(i).normal(i + 1, 0.04, len(vals))
            ax.scatter(x_jitter, vals, color=colors_6[i], s=3.5, alpha=0.6, edgecolors="none", zorder=3)
            med_val = np.median(vals)
            ax.text(i + 1, med_val + 5 if med_val > 0 else 5, f"med={med_val:.1f}",
                    ha="center", va="bottom", fontsize=5.0, color="black", fontweight="bold")
        
    ax.set_ylabel("Biomass Index for megafauna >100 kg ($B_{H,>100}$, relative units)")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5, zorder=1)
    ax.tick_params(axis='x', labelsize=4.8)
    panel_label(ax, "C")

    # ── (D) Individual Body Mass Density ────────────────────────────────────
    ax = fig.add_subplot(gs[2, 0])
    from scipy.stats import gaussian_kde
    
    usable_det = det[
        (det["taxon_quality"].isin(["species", "genus", "family"])) &
        (~det["taxon_quality"].isin(["blank", "human", "domestic"])) &
        (det["body_mass_kg"].notna()) &
        (det["body_mass_kg"] > 0)
    ].copy()
    
    for idx_reg, reg in enumerate(["Congo", "Amazon", "SE_Asia"]):
        sub_det = usable_det[usable_det["region"] == reg]
        log_masses = np.log10(sub_det["body_mass_kg"])
        
        kde = gaussian_kde(log_masses)
        x_vals = np.linspace(log_masses.min() - 0.5, log_masses.max() + 0.5, 200)
        y_vals = kde(x_vals)
        
        ax.plot(x_vals, y_vals, color=PAL[reg], lw=1.2, label=f"{reg} (events={len(sub_det)})")
        ax.fill_between(x_vals, 0, y_vals, color=PAL[reg], alpha=0.15)
        
        med_m = np.median(sub_det["body_mass_kg"])
        ax.axvline(np.log10(med_m), color=PAL[reg], ls="--", lw=0.8, alpha=0.7)
        ax.text(np.log10(med_m), ax.get_ylim()[1] * (0.9 - 0.15 * idx_reg), f" {med_m:.1f} kg",
                color=PAL[reg], fontsize=5.5, fontweight="bold", ha="left" if reg=="Congo" else "right")
        
    ax.set_xlabel("log₁₀ species body mass (kg)")
    ax.set_ylabel("Probability density")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    ax.legend(frameon=False, loc="upper right")
    panel_label(ax, "D")

    # ── (E) Megafauna Biomass Fraction ──────────────────────────────────────
    ax = fig.add_subplot(gs[2, 1])
    congo_pres_frac = congo[congo["B_H_gt50"] > 0]["megafauna_fraction"].values
    congo_abs_frac = congo[congo["B_H_gt50"] == 0]["megafauna_fraction"].values
    amazon_pres_frac = amazon[amazon["B_H_gt50"] > 0]["megafauna_fraction"].values
    amazon_abs_frac = amazon[amazon["B_H_gt50"] == 0]["megafauna_fraction"].values
    seasia_pres_frac = seasia[seasia["B_H_gt50"] > 0]["megafauna_fraction"].values
    seasia_abs_frac = seasia[seasia["B_H_gt50"] == 0]["megafauna_fraction"].values
    
    box_data_frac = [congo_pres_frac, congo_abs_frac, amazon_pres_frac, amazon_abs_frac, seasia_pres_frac, seasia_abs_frac]
    labels_frac = [
        f"Congo\nPres\n(n={len(congo_pres_frac)})",
        f"Congo\nAbs\n(n={len(congo_abs_frac)})",
        f"Amazon\nPres\n(n={len(amazon_pres_frac)})",
        f"Amazon\nAbs\n(n={len(amazon_abs_frac)})",
        f"SE_Asia\nPres\n(n={len(seasia_pres_frac)})",
        f"SE_Asia\nAbs\n(n={len(seasia_abs_frac)})"
    ]
    
    bp_frac = ax.boxplot(box_data_frac, tick_labels=labels_frac, patch_artist=True,
                          widths=0.5, showfliers=False, zorder=2)
    
    for patch, color in zip(bp_frac["boxes"], colors_6):
        patch.set_facecolor(color)
        patch.set_alpha(0.65)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.5)
    for element in ["whiskers", "caps", "medians"]:
        plt.setp(bp_frac[element], color="black", lw=0.5)
        
    for i, vals in enumerate(box_data_frac):
        if len(vals) > 0:
            x_jitter = np.random.default_rng(i).normal(i + 1, 0.04, len(vals))
            ax.scatter(x_jitter, vals, color=colors_6[i], s=3.5, alpha=0.6, edgecolors="none", zorder=3)
            
            mean_val = np.mean(vals)
            ax.text(i + 1, mean_val + 2 if mean_val > 0 else 2, f"mean={mean_val:.1f}%",
                    ha="center", va="bottom", fontsize=5.0, color="black", fontweight="bold")
        
    ax.set_ylabel("Megafaunal Biomass Fraction (%)")
    ax.set_ylim(-2, 105)
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5, zorder=1)
    ax.tick_params(axis='x', labelsize=4.8)
    panel_label(ax, "E")

    # Save
    for ext in ["pdf", "png"]:
        fig.savefig(fig_dir / f"camera_traps_fig03_bodysize.{ext}", dpi=300)
    plt.close(fig)
    print(f"Saved landscape cluster-level camera_traps_fig03_bodysize.pdf/png.")


def make_figure4(cluster_metrics: pd.DataFrame, fig_dir: Path):
    """
    4-panel inspection figure plotting biomass indices and species richness against trap days (cumulative effort)
    to check for sampling effort bias across spatial clusters.
    (A) total standing heterotroph biomass index (B_H_index) vs trap_days (log-log scale)
    (B) standing biomass index of animals > 50 kg (B_H_gt50) vs trap_days (log-linear / symlog scale)
    (C) standing biomass index of megafauna > 100 kg (B_H_gt100) vs trap_days (log-linear / symlog scale)
    (D) species richness (n_species) vs trap_days (log-linear scale)
    """
    from scipy.stats import spearmanr
    
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.85), constrained_layout=True)
    robust_metrics = cluster_metrics[cluster_metrics["trap_days"] > 100.0]

    # ── (A) total standing heterotroph biomass index (B_H_index) vs trap_days ────
    ax = axes[0, 0]
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub = cluster_metrics[cluster_metrics["region"] == region]
        ax.scatter(sub["trap_days"], sub["B_H_index"],
                   c=PAL[region], s=16, alpha=0.7, edgecolors="white", linewidths=0.2,
                   label=region, zorder=3)
        
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cumulative effort (trap-days)")
    ax.set_ylabel("Biomass Index ($B_H$)")
    ax.axvline(100.0, color="#D32F2F", linestyle="--", linewidth=0.8, alpha=0.8, zorder=2, label="Threshold (100 days)")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    
    # Compute spearman correlations (Full vs Filtered)
    r_full, p_full = spearmanr(cluster_metrics["trap_days"], cluster_metrics["B_H_index"])
    r_rob, p_rob = spearmanr(robust_metrics["trap_days"], robust_metrics["B_H_index"])
    
    stat_text = (
        f"Full (n={len(cluster_metrics)}):\n"
        f"  $r_S$ = {r_full:.3f} ($P$ = {p_full:.3g})\n"
        f"Robust (n={len(robust_metrics)}):\n"
        f"  $r_S$ = {r_rob:.3f} ($P$ = {p_rob:.3g})"
    )
    ax.text(0.05, 0.95, stat_text, 
            transform=ax.transAxes, fontsize=4.8, color="black", fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5),
            verticalalignment="top")
    panel_label(ax, "A")
    ax.legend(frameon=False, loc="lower right", fontsize=5.5)

    # ── (B) standing biomass index of animals > 50 kg (B_H_gt50) vs trap_days ────
    ax = axes[0, 1]
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub = cluster_metrics[cluster_metrics["region"] == region]
        ax.scatter(sub["trap_days"], sub["B_H_gt50"],
                   c=PAL[region], s=16, alpha=0.7, edgecolors="white", linewidths=0.2,
                   zorder=3)
        
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("Cumulative effort (trap-days)")
    ax.set_ylabel("Biomass Index for animals >50 kg ($B_{H,>50}$)")
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.axvline(100.0, color="#D32F2F", linestyle="--", linewidth=0.8, alpha=0.8, zorder=2)
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    
    r_full, p_full = spearmanr(cluster_metrics["trap_days"], cluster_metrics["B_H_gt50"])
    r_rob, p_rob = spearmanr(robust_metrics["trap_days"], robust_metrics["B_H_gt50"])
    
    stat_text = (
        f"Full (n={len(cluster_metrics)}):\n"
        f"  $r_S$ = {r_full:.3f} ($P$ = {p_full:.3g})\n"
        f"Robust (n={len(robust_metrics)}):\n"
        f"  $r_S$ = {r_rob:.3f} ($P$ = {p_rob:.3g})"
    )
    ax.text(0.05, 0.95, stat_text, 
            transform=ax.transAxes, fontsize=4.8, color="black", fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5),
            verticalalignment="top")
    panel_label(ax, "B")

    # ── (C) standing biomass index of megafauna > 100 kg (B_H_gt100) vs trap_days ────
    ax = axes[1, 0]
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub = cluster_metrics[cluster_metrics["region"] == region]
        ax.scatter(sub["trap_days"], sub["B_H_gt100"],
                   c=PAL[region], s=16, alpha=0.7, edgecolors="white", linewidths=0.2,
                   zorder=3)
        
    ax.set_xscale("log")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("Cumulative effort (trap-days)")
    ax.set_ylabel("Biomass Index for megafauna >100 kg ($B_{H,>100}$)")
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.axvline(100.0, color="#D32F2F", linestyle="--", linewidth=0.8, alpha=0.8, zorder=2)
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    
    r_full, p_full = spearmanr(cluster_metrics["trap_days"], cluster_metrics["B_H_gt100"])
    r_rob, p_rob = spearmanr(robust_metrics["trap_days"], robust_metrics["B_H_gt100"])
    
    stat_text = (
        f"Full (n={len(cluster_metrics)}):\n"
        f"  $r_S$ = {r_full:.3f} ($P$ = {p_full:.3g})\n"
        f"Robust (n={len(robust_metrics)}):\n"
        f"  $r_S$ = {r_rob:.3f} ($P$ = {p_rob:.3g})"
    )
    ax.text(0.05, 0.95, stat_text, 
            transform=ax.transAxes, fontsize=4.8, color="black", fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5),
            verticalalignment="top")
    panel_label(ax, "C")

    # ── (D) species richness (n_species) vs trap_days ────
    ax = axes[1, 1]
    for region in ["Congo", "Amazon", "SE_Asia"]:
        sub = cluster_metrics[cluster_metrics["region"] == region]
        ax.scatter(sub["trap_days"], sub["n_species"],
                   c=PAL[region], s=16, alpha=0.7, edgecolors="white", linewidths=0.2,
                   zorder=3)
        
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Cumulative effort (trap-days)")
    ax.set_ylabel("Cluster taxon richness ($S$)")
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.axvline(100.0, color="#D32F2F", linestyle="--", linewidth=0.8, alpha=0.8, zorder=2)
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    
    r_full, p_full = spearmanr(cluster_metrics["trap_days"], cluster_metrics["n_species"])
    r_rob, p_rob = spearmanr(robust_metrics["trap_days"], robust_metrics["n_species"])
    
    stat_text = (
        f"Full (n={len(cluster_metrics)}):\n"
        f"  $r_S$ = {r_full:.3f} ($P$ = {p_full:.3g})\n"
        f"Robust (n={len(robust_metrics)}):\n"
        f"  $r_S$ = {r_rob:.3f} ($P$ = {p_rob:.3g})"
    )
    ax.text(0.05, 0.95, stat_text, 
            transform=ax.transAxes, fontsize=4.8, color="black", fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=1.5),
            verticalalignment="top")
    panel_label(ax, "D")

    # Save
    for ext in ["pdf", "png"]:
        fig.savefig(fig_dir / f"camera_traps_fig04_bias_check.{ext}", dpi=300)
    plt.close(fig)
    print(f"Saved landscape cluster-level camera_traps_fig04_bias_check.pdf/png.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate PNAS-style multi-panel figures from pre-processed camera trap datasets."
    )
    parser.add_argument("--fig-dir", type=Path, default=FIG_DIR)
    args = parser.parse_args()

    args.fig_dir.mkdir(parents=True, exist_ok=True)
    set_pnas_style()

    # Load pre-processed datasets
    robust_det_path = OUTPUT_DIR / "camera_traps_robust_detections.csv"
    robust_metrics_path = OUTPUT_DIR / "camera_traps_cluster_level_metrics_robust.csv"
    full_metrics_path = OUTPUT_DIR / "camera_traps_cluster_level_metrics.csv"
    
    if not robust_det_path.exists() or not robust_metrics_path.exists() or not full_metrics_path.exists():
        raise FileNotFoundError(
            "Pre-processed outputs are missing. Please run code/process_camera_traps.py first."
        )

    print("Loading pre-processed camera trapping data...")
    robust_det = pd.read_csv(robust_det_path)
    robust_cluster_metrics = pd.read_csv(robust_metrics_path)
    cluster_metrics = pd.read_csv(full_metrics_path)

    print("Loading pre-computed spatial cluster geometries from GeoJSON...")
    geojson_path = OUTPUT_DIR / "camera_traps_robust_buffered_mcps.geojson"
    if not geojson_path.exists():
        raise FileNotFoundError(
            f"Pre-processed GeoJSON is missing: {geojson_path}. Please run code/process_camera_traps.py first."
        )
    with open(geojson_path, "w" if False else "r") as f:
        geojson_data = json.load(f)
    
    geojson_features = []
    for feat in geojson_data["features"]:
        geojson_features.append({
            "cluster_id": feat["properties"]["cluster_id"],
            "region": feat["properties"]["region"],
            "geometry": shape(feat["geometry"]),
            "trap_days": feat["properties"]["trap_days"]
        })

    # Load config file for min_trap_days
    config_path = Path(__file__).resolve().parent / "config.json"
    if config_path.exists():
        with open(config_path, "r") as f:
            config = json.load(f)
        min_trap_days = config.get("clustering", {}).get("min_trap_days", 10)
    else:
        min_trap_days = 10

    print(f"\nTotal GEDI-valid clusters: {len(cluster_metrics)} | Robust clusters (>= {min_trap_days} trap-days): {len(robust_cluster_metrics)}")

    print("\nGenerating Figure 1: Detections & Diversity Comparison (with Buffered MCPs)...")
    make_figure1(robust_det, robust_cluster_metrics, geojson_features, args.fig_dir)

    print("\nGenerating Figure 2: Biophysical Scaling & Energetics Comparison...")
    make_figure2(robust_cluster_metrics, args.fig_dir)

    print("\nGenerating Figure 3: Vertebrate Body Size & Megafauna Comparison...")
    make_figure3(robust_det, robust_cluster_metrics, geojson_features, args.fig_dir)

    print("\nGenerating Figure 4: Effort-Bias Diagnostics Comparison...")
    make_figure4(cluster_metrics, args.fig_dir)

    print(f"\nAll cluster-level PNAS camera trap figures successfully generated and saved to {args.fig_dir}/")


if __name__ == "__main__":
    main()
