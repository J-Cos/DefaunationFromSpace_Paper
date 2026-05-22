#!/usr/bin/env python3
"""
visualise_camera_traps.py
=========================
Generates publication-quality (PNAS-style) multi-panel figures comparing
vertebrate community structure and biophysical scaling across the
Amazon and Congo Basins.

Figure 1: Camera Trap Detections & Diversity Comparison
    (A) Maps of camera trap deployments in Congo (left) and Amazon (right)
    (B) Rank-abundance curves of species in both basins
    (C) Taxonomic composition by order (proportional horizontal bars)
    (D) Deployment-level RAI distributions (box/violin plot)

Figure 2: Biophysical Scaling & Energetics Comparison
    (A) log B_H vs log M_H scatter plot with fitted scaling lines
    (B) Species richness vs log B_H scatter plot
    (C) Standing heterotroph biomass (B_H) index distribution
    (D) Standing heterotroph metabolism (M_H) index distribution

Output: - figures/camera_traps_fig01_detections.pdf/png
        - figures/camera_traps_fig02_biophysical.pdf/png
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from pathlib import Path

# ── Paths ───────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
OUTPUT_DIR = PROJECT_DIR / "outputs"
FIG_DIR = PROJECT_DIR / "figures"
FIG_DIR.mkdir(exist_ok=True)

# ── PNAS style configuration ───────────────────────────────────────────────
# PNAS: single column = 3.42 in (8.7 cm), double column = 7.0 in (17.8 cm)
SINGLE_COL = 3.42
DOUBLE_COL = 7.0

# Colors - muted, professional, colorblind-safe
PAL = {
    "Congo":   "#2166AC",  # Blue for Congo
    "Amazon":  "#B2182B",  # Red for Amazon
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


# ── Figure 1: Camera Trap Detections & Diversity ───────────────────────────

def make_figure1(det: pd.DataFrame, fig_dir: Path):
    """
    4-panel summary of camera trap detections and diversity across basins.
    (A) Deployment maps (Congo & Amazon)
    (B) Rank-abundance comparison
    (C) Proportional taxonomic order composition
    (D) RAI distribution across deployments
    """
    # ── Map Ingestion & Import Cartopy ──────────────────────────────────────
    try:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        cartopy_available = True
    except ImportError:
        cartopy_available = False

    # Unique deployments
    dep = (det.groupby(["region", "project_name", "deployment_id"])
           .agg(lon=("longitude", "first"),
                lat=("latitude", "first"),
                trap_days=("trap_days", "first"),
                n_detections=("n_detections", "sum"))
           .reset_index())

    fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.95))
    
    # ── Gridspec Layout ─────────────────────────────────────────────────────
    # Top half: Map A (split into two subpanels for Congo and Amazon)
    # Bottom half: Panels B, C, D
    gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 1.1], hspace=0.35, wspace=0.35)

    # ── (A) Geographic Maps (Left = Congo, Right = Amazon) ──────────────────
    # We will create a nested subgridspec for the top row to show two maps side-by-side
    gs_maps = gs[0, :].subgridspec(1, 2, wspace=0.15)
    
    regions_info = [
        {"name": "Congo", "extent": [9.0, 17.5, -2.5, 5.5], "proj_lon": 13.0, "proj_lat": 1.5, "gs": gs_maps[0], "lbl": "A1: Congo deployments"},
        {"name": "Amazon", "extent": [-63.5, -45.5, -8.5, 3.5], "proj_lon": -54.5, "proj_lat": -2.5, "gs": gs_maps[1], "lbl": "A2: Amazon deployments"}
    ]

    for idx, reg in enumerate(regions_info):
        sub_dep = dep[dep["region"] == reg["name"]]
        
        if cartopy_available:
            ax_map = fig.add_subplot(reg["gs"], projection=ccrs.PlateCarree())
            ax_map.set_extent(reg["extent"], crs=ccrs.PlateCarree())
            
            # Basemap features
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

        # Plot deployments sized by trap-days
        sizes = np.clip(sub_dep["trap_days"] / 6.0, 4, 80)
        color = PAL[reg["name"]]
        
        if cartopy_available:
            sc = ax_map.scatter(
                sub_dep["lon"], sub_dep["lat"],
                s=sizes, c=color, alpha=0.6,
                edgecolors="white", linewidths=0.2,
                transform=ccrs.PlateCarree(), zorder=3
            )
        else:
            sc = ax_map.scatter(
                sub_dep["lon"], sub_dep["lat"],
                s=sizes, c=color, alpha=0.6,
                edgecolors="white", linewidths=0.2, zorder=3
            )
            
        ax_map.set_title(f"{reg['name']} Basin (n={len(sub_dep)})", fontsize=7.5, pad=3)
        if idx == 0:
            panel_label(ax_map, "A", x=-0.08, y=1.05)
            
            # Map scale/sizes legend
            for td, sz in [(30, 5), (150, 25), (365, 60)]:
                ax_map.scatter([], [], s=sz, c=PAL["grey"], alpha=0.6,
                               edgecolors="white", linewidths=0.2,
                               label=f"{td} d")
            ax_map.legend(fontsize=4.5, frameon=True, fancybox=False,
                          framealpha=0.9, edgecolor="#E0E0E0",
                          loc="lower left", title="Effort (trap-days)", title_fontsize=4.5)

    # ── (B) Rank-Abundance Curves ───────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 0])
    
    for region in ["Congo", "Amazon"]:
        sub_det = det[(det["region"] == region) & (det["taxon_quality"] == "species")].copy()
        
        # Aggregate total detections by species binomial
        sp_counts = (sub_det.groupby(["genus", "species", "common_name"])["n_detections"]
                     .sum().sort_values(ascending=False).reset_index())
        sp_counts["rank"] = range(1, len(sp_counts) + 1)
        
        ax.plot(sp_counts["rank"], sp_counts["n_detections"], 
                color=PAL[region], lw=1.2, marker="o", ms=2.5, alpha=0.8,
                label=f"{region} (S={len(sp_counts)})")
        
        # Annotate top species
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
    for region in ["Congo", "Amazon"]:
        sub_det = det[(det["region"] == region) & (det["taxon_quality"].isin(["species", "genus", "family"]))].copy()
        
        # Exclude blanks, humans, domestic
        sub_det = sub_det[~sub_det["taxon_quality"].isin(["blank", "human", "domestic"])]
        
        # Sum detections by order
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
    
    # Order categories logically so legend matches colors nicely
    legend_orders = [o for o in ORDER_COLOURS if o in df_pivot.columns]
    if "Other" in df_pivot.columns and "Other" not in legend_orders:
        legend_orders.append("Other")
    df_pivot = df_pivot[legend_orders]

    left = np.zeros(2)
    regions_idx = ["Congo", "Amazon"]
    
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

    # ── (D) RAI Distribution ────────────────────────────────────────────────
    ax = fig.add_subplot(gs[1, 2])
    
    # Calculate RAI per deployment (sum of wild species only)
    wild_det = det[
        (det["taxon_quality"].isin(["species", "genus", "family"])) &
        (~det["taxon_quality"].isin(["blank", "human", "domestic"]))
    ].copy()
    
    site_rai = wild_det.groupby(["region", "deployment_id"])["RAI"].sum().reset_index()
    
    # Create simple boxplots
    box_data = [
        site_rai[site_rai["region"] == "Congo"]["RAI"].values,
        site_rai[site_rai["region"] == "Amazon"]["RAI"].values
    ]
    
    bp = ax.boxplot(box_data, labels=["Congo", "Amazon"], patch_artist=True,
                    widths=0.45, showfliers=False, zorder=2)
    
    colors = [PAL["Congo"], PAL["Amazon"]]
    for patch, color in zip(bp["boxes"], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.6)
        patch.set_edgecolor("black")
        patch.set_linewidth(0.6)
        
    for element in ["whiskers", "caps", "medians"]:
        plt.setp(bp[element], color="black", lw=0.6)
        
    # Jitter raw points
    for i, region in enumerate(["Congo", "Amazon"]):
        vals = box_data[i]
        x_jitter = np.random.default_rng(i).normal(i + 1, 0.04, len(vals))
        ax.scatter(x_jitter, vals, color=colors[i], s=5, alpha=0.35, edgecolors="none", zorder=3)
        
        # Add median label
        median_val = np.median(vals)
        ax.text(i + 1, median_val * 1.15 if median_val > 0 else 5, f"med={median_val:.1f}",
                ha="center", va="bottom", fontsize=5, color="black", fontweight="bold")
        
    ax.set_ylabel("Wild heterotroph RAI\n(detections / 100 trap-days)")
    ax.set_yscale("symlog", linthresh=1.0)
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5, zorder=1)
    panel_label(ax, "D")

    # ── Save Figures ────────────────────────────────────────────────────────
    for ext in ["pdf", "png"]:
        fig.savefig(fig_dir / f"camera_traps_fig01_detections.{ext}", dpi=300)
    plt.close(fig)
    print(f"Saved camera_traps_fig01_detections.pdf/png in PNAS format.")


# ── Figure 2: Biophysical Scaling & Energetics ─────────────────────────────

def make_figure2(sites: pd.DataFrame, fig_dir: Path):
    """
    4-panel biophysical scaling and energetics summary across basins.
    (A) log B_H vs log M_H scatter plot
    (B) Species richness vs log B_H scatter plot
    (C) Standing heterotroph biomass (B_H) index distribution
    (D) Standing heterotroph metabolism (M_H) index distribution
    """
    fig, axes = plt.subplots(2, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.85), constrained_layout=True)

    # Filter usable sites (where biomass and metabolism indices are > 0)
    valid = sites[(sites["B_H_index"] > 0) & (sites["M_H_index"] > 0)].copy()

    # ── (A) log B_H vs log M_H Scatter ──────────────────────────────────────
    ax = axes[0, 0]
    
    # Plot by region
    for region in ["Congo", "Amazon"]:
        sub = valid[valid["region"] == region]
        ax.scatter(np.log10(sub["B_H_index"]), np.log10(sub["M_H_index"]),
                   c=PAL[region], s=12, alpha=0.6, edgecolors="white", linewidths=0.2,
                   label=region, zorder=3)
        
    # Fit pooled line
    logB = np.log10(valid["B_H_index"])
    logM = np.log10(valid["M_H_index"])
    slope, intercept = np.polyfit(logB, logM, 1)
    
    xfit = np.linspace(logB.min(), logB.max(), 50)
    yfit = slope * xfit + intercept
    ax.plot(xfit, yfit, color=PAL["grey"], ls="--", lw=0.9, zorder=2)
    
    # Spearman r
    corr = valid["B_H_index"].corr(valid["M_H_index"], method="spearman")
    ax.text(0.05, 0.92, f"slope = {slope:.3f}\n$r_S$ = {corr:.3f}", 
            transform=ax.transAxes, fontsize=5.5, color="black", fontweight="bold",
            bbox=dict(facecolor="white", alpha=0.8, edgecolor="none", pad=1))
    
    ax.set_xlabel("log₁₀ standing Biomass Index ($B_H$)")
    ax.set_ylabel("log₁₀ standing Metabolism Index ($M_H$)")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    ax.legend(frameon=False, loc="lower right")
    panel_label(ax, "A")

    # ── (B) Species Richness vs log B_H ─────────────────────────────────────
    ax = axes[0, 1]
    
    for region in ["Congo", "Amazon"]:
        sub = valid[valid["region"] == region]
        ax.scatter(sub["n_species"], np.log10(sub["B_H_index"]),
                   c=PAL[region], s=12, alpha=0.6, edgecolors="white", linewidths=0.2,
                   zorder=3)
        
        # Fit trend line per region
        if len(sub) > 2:
            s_slope, s_int = np.polyfit(sub["n_species"], np.log10(sub["B_H_index"]), 1)
            x_s = np.linspace(sub["n_species"].min(), sub["n_species"].max(), 50)
            ax.plot(x_s, s_slope * x_s + s_int, color=PAL[region], ls=":", lw=0.8, zorder=2)

    ax.set_xlabel("Taxon richness ($S$)")
    ax.set_ylabel("log₁₀ standing Biomass Index ($B_H$)")
    ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    panel_label(ax, "B")

    # ── (C) Standing Biomass (B_H) Distribution ────────────────────────────
    ax = axes[1, 0]
    
    logB_congo = np.log10(valid[valid["region"] == "Congo"]["B_H_index"])
    logB_amazon = np.log10(valid[valid["region"] == "Amazon"]["B_H_index"])
    
    # Combined hist with overlapping transparency
    bins = np.histogram(np.hstack((logB_congo, logB_amazon)), bins=20)[1]
    
    ax.hist(logB_congo, bins=bins, color=PAL["Congo"], alpha=0.55, edgecolor=PAL["Congo"],
            linewidth=0.4, label="Congo", density=False)
    ax.hist(logB_amazon, bins=bins, color=PAL["Amazon"], alpha=0.55, edgecolor=PAL["Amazon"],
            linewidth=0.4, label="Amazon", density=False)
    
    # Add median vertical lines
    med_c = np.median(logB_congo)
    med_a = np.median(logB_amazon)
    
    ax.axvline(med_c, color=PAL["Congo"], ls="--", lw=0.8, zorder=4)
    ax.axvline(med_a, color=PAL["Amazon"], ls="--", lw=0.8, zorder=4)
    
    # Label medians (10^med)
    ax.text(med_c + 0.05, ax.get_ylim()[1] * 0.85, f"med={10**med_c:.1f}", 
            color=PAL["Congo"], fontsize=5, fontweight="bold")
    ax.text(med_a - 0.05, ax.get_ylim()[1] * 0.70, f"med={10**med_a:.1f}", 
            color=PAL["Amazon"], fontsize=5, fontweight="bold", ha="right")
    
    ax.set_xlabel("log₁₀ standing Biomass Index ($B_H$)")
    ax.set_ylabel("Number of sites")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    panel_label(ax, "C")

    # ── (D) Standing Metabolism (M_H) Distribution ─────────────────────────
    ax = axes[1, 1]
    
    logM_congo = np.log10(valid[valid["region"] == "Congo"]["M_H_index"])
    logM_amazon = np.log10(valid[valid["region"] == "Amazon"]["M_H_index"])
    
    bins_m = np.histogram(np.hstack((logM_congo, logM_amazon)), bins=20)[1]
    
    ax.hist(logM_congo, bins=bins_m, color=PAL["Congo"], alpha=0.55, edgecolor=PAL["Congo"],
            linewidth=0.4, label="Congo", density=False)
    ax.hist(logM_amazon, bins=bins_m, color=PAL["Amazon"], alpha=0.55, edgecolor=PAL["Amazon"],
            linewidth=0.4, label="Amazon", density=False)
    
    med_mc = np.median(logM_congo)
    med_ma = np.median(logM_amazon)
    
    ax.axvline(med_mc, color=PAL["Congo"], ls="--", lw=0.8, zorder=4)
    ax.axvline(med_ma, color=PAL["Amazon"], ls="--", lw=0.8, zorder=4)
    
    ax.text(med_mc + 0.05, ax.get_ylim()[1] * 0.85, f"med={10**med_mc:.1f}", 
            color=PAL["Congo"], fontsize=5, fontweight="bold")
    ax.text(med_ma - 0.05, ax.get_ylim()[1] * 0.70, f"med={10**med_ma:.1f}", 
            color=PAL["Amazon"], fontsize=5, fontweight="bold", ha="right")
    
    ax.set_xlabel("log₁₀ standing Metabolism Index ($M_H$)")
    ax.set_ylabel("Number of sites")
    ax.grid(True, linestyle="--", linewidth=0.2, color="#E0E0E0", alpha=0.5)
    panel_label(ax, "D")

    # ── Save Figures ────────────────────────────────────────────────────────
    for ext in ["pdf", "png"]:
        fig.savefig(fig_dir / f"camera_traps_fig02_biophysical.{ext}", dpi=300)
    plt.close(fig)
    print(f"Saved camera_traps_fig02_biophysical.pdf/png in PNAS format.")


def main():
    parser = argparse.ArgumentParser(
        description="Generate PNAS-style publication-quality multi-panel figures comparing Amazon and Congo camera trap data."
    )
    parser.add_argument("--detections", type=Path, default=OUTPUT_DIR / "camera_traps_joint_detections.csv")
    parser.add_argument("--metrics", type=Path, default=OUTPUT_DIR / "camera_traps_joint_metrics.csv")
    parser.add_argument("--fig-dir", type=Path, default=FIG_DIR)
    args = parser.parse_args()

    args.fig_dir.mkdir(parents=True, exist_ok=True)
    set_pnas_style()

    print("Loading joint camera trapping datasets...")
    det = pd.read_csv(args.detections)
    metrics = pd.read_csv(args.metrics)

    print("\nGenerating Figure 1: Detections & Diversity Comparison...")
    make_figure1(det, args.fig_dir)

    print("\nGenerating Figure 2: Biophysical Scaling & Energetics Comparison...")
    make_figure2(metrics, args.fig_dir)

    print(f"\nAll camera trap figures successfully generated and saved to {args.fig_dir}/")


if __name__ == "__main__":
    main()
