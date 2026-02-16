import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from itertools import combinations


# ------------------------------------------------------------------
# Basic metrics
# ------------------------------------------------------------------

def compute_clone_table(adata, clone_col, celltype_col):
    return pd.crosstab(
        adata.obs[clone_col],
        adata.obs[celltype_col]
    )


def compute_basic_stats(clone_table):
    cells_per_clone = clone_table.sum(axis=1)

    stats = {
        "total_clones": len(clone_table),
        "total_cells": clone_table.sum().sum(),
        "mean_clone_size": cells_per_clone.mean(),
        "median_clone_size": cells_per_clone.median(),
        "max_clone_size": cells_per_clone.max(),
        "single_cell_clones": (cells_per_clone == 1).sum(),
        "expanded_clones": (cells_per_clone >= 2).sum(),
    }

    stats["cells_per_clone"] = cells_per_clone

    return stats


# ------------------------------------------------------------------
# Diversity metrics
# ------------------------------------------------------------------

def gini_coefficient(x):
    x = np.sort(x)
    n = len(x)
    cumx = np.cumsum(x)
    return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n


def shannon_entropy(series):
    proportions = series / series.sum()
    return -np.sum(proportions * np.log2(proportions + 1e-10))


def simpson_diversity(series):
    proportions = series / series.sum()
    return 1 - np.sum(proportions ** 2)


def compute_diversity_metrics(cells_per_clone):
    gini = gini_coefficient(cells_per_clone.values)
    shannon = shannon_entropy(cells_per_clone)
    max_entropy = np.log2(len(cells_per_clone))
    normalized_entropy = shannon / max_entropy
    simpson = simpson_diversity(cells_per_clone)
    expansion_index = cells_per_clone[cells_per_clone >= 2].sum() / cells_per_clone.sum()

    return {
        "gini": gini,
        "shannon": shannon,
        "normalized_entropy": normalized_entropy,
        "simpson": simpson,
        "expansion_index": expansion_index
    }


def compute_jaccard_matrix(adata, clone_col, celltype_col):
    cell_types = sorted(adata.obs[celltype_col].unique())
    n = len(cell_types)
    matrix = np.zeros((n, n))

    for i in range(n):
        clones_i = set(
            adata[adata.obs[celltype_col] == cell_types[i]].obs[clone_col]
        )
        for j in range(n):
            if i == j:
                matrix[i, j] = 1.0
            elif i < j:
                clones_j = set(
                    adata[adata.obs[celltype_col] == cell_types[j]].obs[clone_col]
                )
                intersection = len(clones_i & clones_j)
                union = len(clones_i | clones_j)
                similarity = intersection / union if union > 0 else 0
                matrix[i, j] = similarity
                matrix[j, i] = similarity

    return cell_types, matrix


def run_clonality_analysis(
    adata,
    clone_col="larry_clone",
    celltype_col="celltype_work"
):
    clone_table = compute_clone_table(adata, clone_col, celltype_col)

    basic = compute_basic_stats(clone_table)
    diversity = compute_diversity_metrics(basic["cells_per_clone"])

    cell_types, jaccard_matrix = compute_jaccard_matrix(
        adata,
        clone_col,
        celltype_col
    )

    results = {
        "basic_stats": basic,
        "diversity": diversity,
        "jaccard_matrix": jaccard_matrix,
        "cell_types": cell_types
    }

    return results



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from itertools import combinations
from scipy import stats as scipy_stats

plt.rcParams.update({
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 20,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16
})



def _compute_clonality_core(adata):
    """
    Compute all clonality statistics needed for Fig.
    """

    # Clone x CellType table
    full_clone_celltype = pd.crosstab(
        adata.obs['larry_clone'],
        adata.obs['CellType_man']
    )

    cells_per_clone = full_clone_celltype.sum(axis=1)
    clone_sizes = cells_per_clone.value_counts().sort_index()

    total_clones = len(full_clone_celltype)
    total_cells = full_clone_celltype.sum().sum()

    clones_with_1_cell = sum(cells_per_clone == 1)
    clones_with_2plus_cells = sum(cells_per_clone >= 2)

    # ---- Gini ----
    def gini_coefficient(x):
        x = np.sort(x)
        n = len(x)
        cumx = np.cumsum(x)
        return (n + 1 - 2 * np.sum(cumx) / cumx[-1]) / n

    gini = gini_coefficient(cells_per_clone.values)

    # ---- Shannon ----
    def shannon_entropy(series):
        proportions = series / series.sum()
        return -np.sum(proportions * np.log2(proportions + 1e-10))

    shannon = shannon_entropy(cells_per_clone)
    max_entropy = np.log2(len(cells_per_clone))
    normalized_entropy = shannon / max_entropy

    # ---- Simpson ----
    def simpson_diversity(series):
        proportions = series / series.sum()
        return 1 - np.sum(proportions ** 2)

    simpson = simpson_diversity(cells_per_clone)

    expansion_index = cells_per_clone[cells_per_clone >= 2].sum() / total_cells

    # ---- Cell type statistics ----
    celltype_clone_counts = {}
    for cell_type in adata.obs['CellType_man'].unique():
        type_cells = adata[adata.obs['CellType_man'] == cell_type]
        type_clones = type_cells.obs['larry_clone'].nunique()

        celltype_clone_counts[cell_type] = {
            'n_cells': len(type_cells),
            'n_clones': type_clones,
            'cells_per_clone': len(type_cells) / type_clones if type_clones > 0 else 0
        }

    # ---- Jaccard ----
    def calculate_jaccard(celltype1, celltype2):
        clones1 = set(adata[adata.obs['CellType_man'] == celltype1].obs['larry_clone'])
        clones2 = set(adata[adata.obs['CellType_man'] == celltype2].obs['larry_clone'])
        intersection = len(clones1.intersection(clones2))
        union = len(clones1.union(clones2))
        return intersection / union if union > 0 else 0

    cell_types = sorted(adata.obs['CellType_man'].unique())
    n_types = len(cell_types)
    jaccard_matrix = np.zeros((n_types, n_types))

    for i in range(n_types):
        for j in range(n_types):
            if i == j:
                jaccard_matrix[i, j] = 1.0
            elif i < j:
                similarity = calculate_jaccard(cell_types[i], cell_types[j])
                jaccard_matrix[i, j] = similarity
                jaccard_matrix[j, i] = similarity

    return {
        "cells_per_clone": cells_per_clone,
        "clone_sizes": clone_sizes,
        "gini": gini,
        "shannon": shannon,
        "normalized_entropy": normalized_entropy,
        "simpson": simpson,
        "expansion_index": expansion_index,
        "celltype_clone_counts": celltype_clone_counts,
        "jaccard_matrix": jaccard_matrix,
        "cell_types": cell_types,
        "n_types": n_types
    }


# ------------------------------------------------------------------
# FIG1
# ------------------------------------------------------------------

def plot_fig1(adata, save_path=None):
    data = _compute_clonality_core(adata)

    plt.figure(figsize=(6, 5))
    plt.bar(data["clone_sizes"].index,
            data["clone_sizes"].values,
            alpha=0.7,
            color='skyblue')

    plt.xlabel('Cells per clone')
    plt.ylabel('Number of clones')
    plt.title('Clone Size Distribution (Linear)')
    plt.grid(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

    plt.close()


# ------------------------------------------------------------------
# FIG2
# ------------------------------------------------------------------

def plot_fig2(adata, save_path=None):
    data = _compute_clonality_core(adata)

    plt.figure(figsize=(6, 5))

    sorted_sizes = np.sort(data["cells_per_clone"].values)
    yvals = np.arange(1, len(sorted_sizes)+1) / len(sorted_sizes)

    plt.plot(sorted_sizes, yvals, 'b-', linewidth=2)

    plt.xlabel('Clone size (cells)')
    plt.ylabel('Cumulative fraction')
    plt.title('Cumulative Clone Size Distribution')
    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

    plt.close()


# ------------------------------------------------------------------
# FIG3
# ------------------------------------------------------------------

def plot_fig3(adata, save_path=None):
    data = _compute_clonality_core(adata)

    plt.figure(figsize=(6, 5))

    plt.loglog(data["clone_sizes"].index,
               data["clone_sizes"].values,
               'bo', alpha=0.6)

    plt.xlabel('Clone size (log scale)')
    plt.ylabel('Number of clones (log scale)')
    plt.title('Clone Size Distribution (Log-Log)')
    plt.grid(True, alpha=0.3, which='both')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

    plt.close()


# ------------------------------------------------------------------
# FIG4
# ------------------------------------------------------------------

def plot_fig4(adata, save_path=None):
    data = _compute_clonality_core(adata)

    plt.figure(figsize=(6, 5))

    cumulative_clone_sizes = np.cumsum(np.sort(data["cells_per_clone"].values))
    cumulative_clone_fraction = cumulative_clone_sizes / cumulative_clone_sizes[-1]
    perfect_line = np.linspace(0, 1, len(cumulative_clone_fraction))

    plt.plot(perfect_line, perfect_line, 'k--', label='Perfect equality')
    plt.plot(perfect_line, cumulative_clone_fraction, 'r-', linewidth=2, label='Actual distribution')
    plt.fill_between(perfect_line, perfect_line, cumulative_clone_fraction, alpha=0.3)

    plt.xlabel('Cumulative fraction of clones')
    plt.ylabel('Cumulative fraction of cells')
    plt.title(f'Lorenz Curve (Gini = {data["gini"]:.3f})')
    plt.legend()
    plt.grid(True, alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

    plt.close()


# ------------------------------------------------------------------
# FIG5
# ------------------------------------------------------------------

def plot_fig5(adata, save_path=None):
    data = _compute_clonality_core(adata)

    plt.figure(figsize=(8, 6))

    celltypes = list(data["celltype_clone_counts"].keys())
    clone_counts = [stats['n_clones'] for stats in data["celltype_clone_counts"].values()]
    cell_counts = [stats['n_cells'] for stats in data["celltype_clone_counts"].values()]

    x = np.arange(len(celltypes))
    width = 0.35

    plt.bar(x - width/2, cell_counts, width,
            label='Cells', color='lightblue')
    plt.bar(x + width/2, clone_counts, width,
            label='Clones', color='lightcoral')

    plt.xlabel('Cell Type')
    plt.ylabel('Count')
    plt.title('Cells vs Clones per Cell Type')
    plt.xticks(x, celltypes, rotation=45, ha='right')
    plt.legend()
    plt.grid(False)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

    plt.close()


# ------------------------------------------------------------------
# FIG6
# ------------------------------------------------------------------

def plot_fig6(adata, save_path=None):
    data = _compute_clonality_core(adata)

    plt.figure(figsize=(8, 6))

    im = plt.imshow(data["jaccard_matrix"],
                    cmap='Reds',
                    vmin=0,
                    vmax=1)

    plt.xticks(range(data["n_types"]),
               [ct[:8] for ct in data["cell_types"]],
               rotation=45, ha='right')
    plt.yticks(range(data["n_types"]),
               data["cell_types"])

    plt.title('Clone Sharing Between Cell Types')
    plt.colorbar(im, label='Jaccard Similarity')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

    plt.close()


# ------------------------------------------------------------------
# plotting
# ------------------------------------------------------------------

def generate_all_figures(adata, output_dir):
    import os
    os.makedirs(output_dir, exist_ok=True)

    plot_fig1(adata, f"{output_dir}/fig1_clone_size_linear.svg")
    plot_fig2(adata, f"{output_dir}/fig2_cumulative_distribution.svg")
    plot_fig3(adata, f"{output_dir}/fig3_loglog_distribution.svg")
    plot_fig4(adata, f"{output_dir}/fig4_lorenz_curve.svg")
    plot_fig5(adata, f"{output_dir}/fig5_celltype_clone_relationship.svg")
    plot_fig6(adata, f"{output_dir}/fig6_clone_sharing_heatmap.svg")

    
# ------------------------------------------------------------------
# report  
# ------------------------------------------------------------------

def generate_full_clone_report(
    adata,
    output_dir,
    dataset_name=None,
    return_data=True
):
    """
    Generate full clonality report (Fig1–Fig6).

    Parameters
    ----------
    adata : AnnData
    output_dir : str
        Directory to save figures.
    dataset_name : str or None
        Optional prefix for figure titles and filenames.
    return_data : bool
        Whether to return computed clonality statistics.

    Returns
    -------
    dict (optional)
        Core clonality statistics.
    """

    import os
    os.makedirs(output_dir, exist_ok=True)

    # ---- compute once ----
    data = _compute_clonality_core(adata)

    # ---- filename prefix ----
    prefix = f"{dataset_name}_" if dataset_name else ""

    # ---- FIGURES ----
    plot_fig1(adata, f"{output_dir}/{prefix}fig1_clone_size_linear.svg")
    plot_fig2(adata, f"{output_dir}/{prefix}fig2_cumulative_distribution.svg")
    plot_fig3(adata, f"{output_dir}/{prefix}fig3_loglog_distribution.svg")
    plot_fig4(adata, f"{output_dir}/{prefix}fig4_lorenz_curve.svg")
    plot_fig5(adata, f"{output_dir}/{prefix}fig5_celltype_clone_relationship.svg")
    plot_fig6(adata, f"{output_dir}/{prefix}fig6_clone_sharing_heatmap.svg")

    # ---- console summary (paper-style metrics) ----
    print("==================================================")
    if dataset_name:
        print(f"Clonality Report: {dataset_name}")
    else:
        print("Clonality Report")

    print("--------------------------------------------------")
    print(f"Total clones: {len(data['cells_per_clone'])}")
    print(f"Gini coefficient: {data['gini']:.4f}")
    print(f"Shannon entropy: {data['shannon']:.4f}")
    print(f"Normalized entropy: {data['normalized_entropy']:.4f}")
    print(f"Simpson diversity: {data['simpson']:.4f}")
    print(f"Expansion index: {data['expansion_index']:.4f}")
    print("==================================================")

    if return_data:
        return data


# ------------------------------------------------------------------
# Source clone check  
# ------------------------------------------------------------------

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def quick_clone_check(adata, column_to_check, cell_type_A, cell_type_B):
    """Quick one-line results"""
    clones_A = set(adata.obs[adata.obs[column_to_check] == cell_type_A]['larry_clone'].dropna())
    clones_B = set(adata.obs[adata.obs[column_to_check] == cell_type_B]['larry_clone'].dropna())
    shared = clones_A.intersection(clones_B)
    
    print(f"A clones: {len(clones_A)} | B clones: {len(clones_B)} | Shared: {len(shared)} | Prop A goes to B: {len(shared)/len(clones_A) if clones_A else 0:.1%}")
    return len(clones_A), len(clones_B), len(shared)


def analyze_multiple_to_B(adata, column_to_check, source_cell_types, target_B):

    # Get clones for target B
    clones_B = set(adata.obs[adata.obs[column_to_check] == target_B]['larry_clone'].dropna())
    n_clones_B = len(clones_B)
    
    print(f"ANALYSIS: How many clones of different sources go to {target_B}?")
    print("=" * 70)
    print(f"Target cell type: {target_B} ({n_clones_B} clones)")
    print("Metric: (shared clones between source and B) / (source clones)")
    print("=" * 70)
    
    results = []
    
    for source_A in source_cell_types:
        # Get clones for source A
        clones_A = set(adata.obs[adata.obs[column_to_check] == source_A]['larry_clone'].dropna())
        n_clones_A = len(clones_A)
        
        # Find shared clones
        shared_clones = clones_A.intersection(clones_B)
        n_shared = len(shared_clones)
        
        # Calculate percentage: shared/source_A
        percentage = (n_shared / n_clones_A * 100) if n_clones_A > 0 else 0
        
        # Get clones that actually contain both types
        clones_with_both = []
        for clone in shared_clones:
            clone_cell_types = adata.obs[adata.obs['larry_clone'] == clone][column_to_check]
            if (source_A in clone_cell_types.values) and (target_B in clone_cell_types.values):
                clones_with_both.append(clone)
        
        results.append({
            'Source': source_A,
            'Target': target_B,
            'Source Clones': n_clones_A,
            'Target Clones': n_clones_B,
            'Shared Clones': n_shared,
            'Clones with Both': len(clones_with_both),
            'Percentage (shared/source)': percentage,
            'Ratio (shared/source)': f"{n_shared}/{n_clones_A}",
            'Shared Clone IDs': list(shared_clones),
            'Clones with Both IDs': clones_with_both
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results)
    
    # Sort by Percentage (decreasing order) - MAIN SORTING
    results_df = results_df.sort_values('Percentage (shared/source)', ascending=False)
    
    # Print results table
    print("\nRESULTS (sorted by percentage, highest first):")
    print("=" * 80)
    
    display_df = results_df[['Source', 'Source Clones', 'Shared Clones', 
                           'Percentage (shared/source)', 'Ratio (shared/source)']].copy()
    display_df['Percentage (shared/source)'] = display_df['Percentage (shared/source)'].round(2)
    
    print(display_df.to_string(index=False))


def visualize_sources_to_B(adata, column_to_check, source_cell_types, target_B, fig_title, figsize=(10, 6)):
    
    # Get analysis results
    results_df = analyze_multiple_to_B(adata, column_to_check, source_cell_types, target_B)
    
    # Create figure with just one plot
    fig, ax = plt.subplots(figsize=figsize)
    
    # Get data for plotting
    sources = results_df['Source']
    percentages = results_df['Percentage (shared/source)']
    
    # Create horizontal bar chart
    y_pos = np.arange(len(sources))
    bars = ax.barh(y_pos, percentages, color='lightblue', edgecolor='navy', alpha=0.8)
    
    # Customize y-axis
    ax.set_yticks(y_pos)
    ax.set_yticklabels(sources)
    ax.invert_yaxis()  # Highest percentage at the top
    
    # Customize x-axis
    ax.set_xlabel(f'% of source clones that contain {target_B}', fontsize=12)
    ax.set_xlim([0, max(percentages) * 1.15])  # Add some padding
    
    # Add title
    ax.set_title(fig_title, # : Percentage of source clones shared
                fontsize=14, pad=20) # fontweight='bold',
    
    # Add percentage labels on bars
    for bar, percentage, source_clones, shared_clones in zip(bars, percentages, 
                                                           results_df['Source Clones'], 
                                                           results_df['Shared Clones']):
        width = bar.get_width()
        
        # Choose label position based on bar width
        if width > 5:  # Enough space inside bar
            label_x = width - 1
            ha = 'right'
            color = 'white'
            fontweight = 'bold'
        else:  # Outside bar
            label_x = width + 0.5
            ha = 'left'
            color = 'black'
            fontweight = 'normal'
        
        # Add percentage label
        ax.text(label_x, bar.get_y() + bar.get_height()/2, 
               f'{width:.1f}%', ha=ha, va='center', 
               color=color, fontweight=fontweight, fontsize=10)
        
        # # Add clone count annotation on the left
        # ax.text(-max(percentages)*0.05, bar.get_y() + bar.get_height()/2,
        #        f"{shared_clones}/{source_clones}", ha='right', va='center',
        #        fontsize=9, color='gray')
    
    # Add grid for better readability
    ax.grid(True, axis='x', alpha=0.3, linestyle='--')
    
    # Add reference lines
    max_perc = max(percentages)
    for perc in [25, 50, 75]:
        if perc <= max_perc:
            ax.axvline(x=perc, color='red', linestyle=':', alpha=0.3, linewidth=1)
    
    # Adjust layout
    plt.tight_layout()
    plt.show()
    
    return results_df


def plot_shared_clones_AB_only(adata, celltype_label, celltype_A, celltype_B, 
                              colorA, colorB, colorC, save_path=None):

    import scanpy as sc
    import matplotlib.pyplot as plt
    import pandas as pd
    import numpy as np
    
    # Recreate clone_type from larry_clone - assign only to relevant cells
    print("Recovering clone_type from larry_clone...")
    adata.obs['clone_type'] = 'Other'
    
    # Identify mixed clones (clones containing both A and B)
    mixed_clone_ids = []
    for clone_id in adata.obs['larry_clone'].dropna().unique():
        clone_mask = adata.obs['larry_clone'] == clone_id
        celltypes_in_clone = adata.obs.loc[clone_mask, celltype_label].unique()
        
        if (celltype_A in celltypes_in_clone) and (celltype_B in celltypes_in_clone):
            mixed_clone_ids.append(clone_id)
    
    # Assign clone types only to relevant cells
    for clone_id in adata.obs['larry_clone'].dropna().unique():
        clone_mask = adata.obs['larry_clone'] == clone_id
        
        if clone_id in mixed_clone_ids:
            # Only assign "Shared" to celltype_A and celltype_B cells in mixed clones
            type_A_mask = clone_mask & (adata.obs[celltype_label] == celltype_A)
            type_B_mask = clone_mask & (adata.obs[celltype_label] == celltype_B)
            adata.obs.loc[type_A_mask, 'clone_type'] = 'Shared'
            adata.obs.loc[type_B_mask, 'clone_type'] = 'Shared'
        else:
            # For non-mixed clones
            celltypes_in_clone = adata.obs.loc[clone_mask, celltype_label].unique()
            
            if len(celltypes_in_clone) == 1:
                if celltypes_in_clone[0] == celltype_A:
                    # Only assign to celltype_A cells
                    type_A_mask = clone_mask & (adata.obs[celltype_label] == celltype_A)
                    adata.obs.loc[type_A_mask, 'clone_type'] = f'{celltype_A}-only'
                elif celltypes_in_clone[0] == celltype_B:
                    # Only assign to celltype_B cells
                    type_B_mask = clone_mask & (adata.obs[celltype_label] == celltype_B)
                    adata.obs.loc[type_B_mask, 'clone_type'] = f'{celltype_B}-only'
    
    # Create the plots
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    
    # Get UMAP coordinates for ALL cells
    all_umap_x = adata.obsm['X_umap'][:, 0]
    all_umap_y = adata.obsm['X_umap'][:, 1]
    
    # Plot 1: All cells as grey background
    axes[0].scatter(all_umap_x, all_umap_y, s=3, alpha=0.1, color='grey', label='All cells')
    
    # Overlay celltype_A cells
    mask_A = adata.obs[celltype_label] == celltype_A
    umap_x_A = adata.obsm['X_umap'][mask_A.values, 0]
    umap_y_A = adata.obsm['X_umap'][mask_A.values, 1]
    axes[0].scatter(umap_x_A, umap_y_A, s=15, alpha=0.7, color=colorA, label=celltype_A)
    
    # Overlay celltype_B cells
    mask_B = adata.obs[celltype_label] == celltype_B
    umap_x_B = adata.obsm['X_umap'][mask_B.values, 0]
    umap_y_B = adata.obsm['X_umap'][mask_B.values, 1]
    axes[0].scatter(umap_x_B, umap_y_B, s=15, alpha=0.7, color=colorB, label=celltype_B)
    
    axes[0].set_title(f'Cell Type ({celltype_A} & {celltype_B} highlighted)')
    axes[0].legend()
    
    # Remove grid lines, frame, and axis labels
    axes[0].grid(False)
    axes[0].set_frame_on(False)
    axes[0].set_xticks([])
    axes[0].set_yticks([])
    axes[0].set_xlabel('')
    axes[0].set_ylabel('')
    
    # Plot 2: All cells as grey background
    axes[1].scatter(all_umap_x, all_umap_y, s=3, alpha=0.1, color='grey', label='All cells')
    
    # Get indices for different clone types - ONLY for celltype_A and celltype_B
    mask_A_only = (adata.obs['clone_type'] == f'{celltype_A}-only')
    mask_B_only = (adata.obs['clone_type'] == f'{celltype_B}-only')
    mask_shared = (adata.obs['clone_type'] == 'Shared')
    
    # Plot celltype_A-only clone cells
    if mask_A_only.sum() > 0:
        umap_x_A_only = adata.obsm['X_umap'][mask_A_only.values, 0]
        umap_y_A_only = adata.obsm['X_umap'][mask_A_only.values, 1]
        axes[1].scatter(umap_x_A_only, umap_y_A_only, s=15, alpha=0.7, 
                       color=colorA, label=f'{celltype_A}-only clones')
    
    # Plot celltype_B-only clone cells
    if mask_B_only.sum() > 0:
        umap_x_B_only = adata.obsm['X_umap'][mask_B_only.values, 0]
        umap_y_B_only = adata.obsm['X_umap'][mask_B_only.values, 1]
        axes[1].scatter(umap_x_B_only, umap_y_B_only, s=15, alpha=0.7, 
                       color=colorB, label=f'{celltype_B}-only clones')
    
    # Plot shared clone cells
    if mask_shared.sum() > 0:
        umap_x_shared = adata.obsm['X_umap'][mask_shared.values, 0]
        umap_y_shared = adata.obsm['X_umap'][mask_shared.values, 1]
        axes[1].scatter(umap_x_shared, umap_y_shared, s=20, alpha=0.8, 
                       color=colorC, edgecolor='black', linewidth=0.5, 
                       label='Shared clones')
    
    axes[1].set_title('Clone Sharing Pattern')
    
    # Move legend to the right (outside the plot)
    axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
    
    # Remove grid lines, frame, and axis labels
    axes[1].grid(False)
    axes[1].set_frame_on(False)
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    axes[1].set_xlabel('')
    axes[1].set_ylabel('')
    
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')
        print(f"Figure saved to {save_path}")
    
    plt.show()