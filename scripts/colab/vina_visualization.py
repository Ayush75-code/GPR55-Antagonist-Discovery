"""
AutoDock Vina Visualization
===========================
Generates publication-quality plots for docking results.
Run after analysis is complete.

Usage: python vina_visualization.py --dir /path/to/results
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import argparse


# ============================================================================
# Configuration
# ============================================================================

POCKET_NAMES = {
    'P0': 'Orthosteric',
    'P1': 'Side_Pocket',
    'P2': 'Interface',
    'P3': 'Allosteric_PPI',
    'P4': 'Lower_Pocket',
    'P5': 'Surface_Groove'
}

LIGAND_DISPLAY_NAMES = {
    'AM251_control': 'AM251_Control',
    'AM251_pubchem': 'AM251_PubChem'
}

POCKETS_ORDER = ['P0', 'P1', 'P2', 'P3', 'P4', 'P5']


def format_number(num):
    return f"{num:,}"


def load_data(results_dir):
    """Load aggregated summary files"""
    aggregated_dir = os.path.join(results_dir, 'aggregated')
    summary_files = glob.glob(os.path.join(aggregated_dir, '*_summary.csv'))
    
    if not summary_files:
        print("⚠️  No summary files found.")
        return None, None
    
    all_results = []
    
    for f in summary_files:
        df = pd.read_csv(f)
        basename = os.path.basename(f).replace('_summary.csv', '')
        parts = basename.split('_P')
        ligand = parts[0]
        pocket_info = 'P' + parts[1] if len(parts) > 1 else ''
        pocket_parts = pocket_info.split('_', 1)
        pocket_id = pocket_parts[0]
        pocket_name = pocket_parts[1] if len(pocket_parts) > 1 else ''
        
        df['Ligand'] = ligand
        df['Pocket'] = pocket_id
        df['Pocket_Name'] = pocket_name
        all_results.append(df)
    
    combined = pd.concat(all_results, ignore_index=True)
    
    # Calculate statistics
    stats = combined.groupby(['Ligand', 'Pocket', 'Pocket_Name'])['Affinity_kcal_mol'].agg([
        ('Best_Affinity', 'min'),
        ('Mean_Affinity', 'mean'),
        ('Median_Affinity', 'median'),
        ('Std_Affinity', 'std'),
        ('Q1', lambda x: x.quantile(0.25)),
        ('Q3', lambda x: x.quantile(0.75)),
        ('Total_Poses', 'count')
    ]).reset_index()
    
    print(f"✓ Loaded {format_number(len(combined))} binding poses")
    
    return stats, combined


def plot_boxplots(combined, output_dir):
    """Generate box plots comparing binding affinities"""
    print("\nGenerating box plots...")
    
    sns.set_style("whitegrid")
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    
    for idx, ligand in enumerate(combined['Ligand'].unique()):
        ax = axes[idx]
        ligand_data = combined[combined['Ligand'] == ligand].copy()
        ligand_data['Pocket'] = pd.Categorical(
            ligand_data['Pocket'], 
            categories=POCKETS_ORDER, 
            ordered=True
        )
        ligand_data = ligand_data.sort_values('Pocket')
        
        box_parts = ax.boxplot(
            [ligand_data[ligand_data['Pocket'] == p]['Affinity_kcal_mol'].values
             for p in POCKETS_ORDER],
            labels=POCKETS_ORDER,
            patch_artist=True,
            showmeans=True
        )
        
        # Color P0 differently
        for i, box in enumerate(box_parts['boxes']):
            if i == 0:  # P0
                box.set_facecolor('lightgreen')
                box.set_alpha(0.7)
            else:
                box.set_facecolor('lightblue')
                box.set_alpha(0.5)
        
        display_name = LIGAND_DISPLAY_NAMES.get(ligand, ligand)
        ax.set_xlabel('Pocket', fontsize=12, fontweight='bold')
        ax.set_ylabel('Binding Affinity (kcal/mol)', fontsize=12, fontweight='bold')
        ax.set_title(f'{display_name}\n(Green = P0 Orthosteric)', fontsize=13, fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'plot1_binding_comparison.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_heatmap(stats, output_dir):
    """Generate publication-quality heatmap"""
    print("Generating heatmap...")
    
    # Map display names
    stats = stats.copy()
    stats['Pocket_Display'] = stats['Pocket'].map(POCKET_NAMES)
    stats['Ligand_Display'] = stats['Ligand'].map(LIGAND_DISPLAY_NAMES)
    stats['Pocket_Display'].fillna(stats['Pocket'], inplace=True)
    stats['Ligand_Display'].fillna(stats['Ligand'], inplace=True)
    
    # Create pivot table
    heatmap_data = stats.pivot(
        index='Ligand_Display',
        columns='Pocket_Display',
        values='Best_Affinity'
    )
    
    # Reorder columns
    column_order = [name for name in POCKET_NAMES.values() if name in heatmap_data.columns]
    heatmap_data = heatmap_data[column_order]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Color map: green (strong) to red (weak)
    cmap = sns.diverging_palette(145, 10, s=85, l=45, as_cmap=True)
    
    sns.heatmap(
        heatmap_data,
        annot=True,
        fmt='.2f',
        cmap=cmap,
        center=-7.0,
        vmin=-10,
        vmax=-4,
        cbar_kws={
            'label': 'Binding Affinity (kcal/mol)',
            'shrink': 0.8,
            'aspect': 30
        },
        linewidths=3,
        linecolor='white',
        square=True,
        ax=ax,
        annot_kws={'fontsize': 13, 'fontweight': 'bold'}
    )
    
    ax.set_xlabel('Target Site', fontsize=15, fontweight='bold', labelpad=12)
    ax.set_ylabel('Ligand', fontsize=15, fontweight='bold', labelpad=12)
    ax.set_title(
        'Best Binding Affinities: Ligand-Target Matrix',
        fontsize=18,
        fontweight='bold',
        pad=25
    )
    
    plt.xticks(rotation=45, ha='right', fontsize=12, fontweight='bold')
    plt.yticks(rotation=0, fontsize=12, fontweight='bold')
    
    # Add caption
    fig.text(
        0.5, 0.02,
        'Lower (more negative) values indicate stronger binding.',
        ha='center',
        fontsize=11,
        style='italic'
    )
    
    plt.tight_layout(rect=[0, 0.08, 1, 0.96])
    output_file = os.path.join(output_dir, 'binding_affinity_matrix.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor='white')
    print(f"✓ Saved: {output_file}")
    plt.close()


def plot_bar_chart(stats, output_dir):
    """Generate bar chart comparison"""
    print("Generating bar chart...")
    
    fig, ax = plt.subplots(figsize=(14, 8))
    x = np.arange(len(POCKETS_ORDER))
    width = 0.35
    
    ligands = stats['Ligand'].unique()
    
    for idx, ligand in enumerate(ligands):
        ligand_data = stats[stats['Ligand'] == ligand].set_index('Pocket').reindex(POCKETS_ORDER)
        display_name = LIGAND_DISPLAY_NAMES.get(ligand, ligand)
        
        bars = ax.bar(
            x + idx * width, 
            ligand_data['Best_Affinity'], 
            width, 
            label=display_name, 
            alpha=0.8
        )
        
        # Highlight P0
        bars[0].set_edgecolor('green')
        bars[0].set_linewidth(3)
    
    ax.set_xlabel('Pocket', fontsize=12, fontweight='bold')
    ax.set_ylabel('Best Binding Affinity (kcal/mol)', fontsize=12, fontweight='bold')
    ax.set_title(
        'Best Binding Affinity by Pocket\n(Green border = P0 Expected Best)', 
        fontsize=14, 
        fontweight='bold'
    )
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(POCKETS_ORDER)
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    output_file = os.path.join(output_dir, 'plot3_validation_bars.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Saved: {output_file}")
    plt.close()


def classify_binding(affinity):
    """Classify binding strength"""
    if affinity < -8.0:
        return 'Strong'
    elif affinity < -6.0:
        return 'Moderate'
    elif affinity < -4.0:
        return 'Weak'
    else:
        return 'Very Weak'


def save_detailed_table(stats, output_dir):
    """Save detailed comparison table"""
    print("Saving detailed comparison table...")
    
    stats = stats.copy()
    stats['Binding_Strength'] = stats['Best_Affinity'].apply(classify_binding)
    
    comparison_table = stats[[
        'Ligand', 'Pocket', 'Pocket_Name',
        'Best_Affinity', 'Mean_Affinity', 'Median_Affinity',
        'Std_Affinity', 'Binding_Strength', 'Total_Poses'
    ]].copy()
    
    comparison_table = comparison_table.sort_values(['Ligand', 'Best_Affinity'])
    
    output_file = os.path.join(output_dir, 'detailed_comparison_table.csv')
    comparison_table.to_csv(output_file, index=False, float_format='%.2f')
    print(f"✓ Saved: {output_file}")
    
    return comparison_table


def main():
    parser = argparse.ArgumentParser(description='Visualize AutoDock Vina results')
    parser.add_argument('--dir', default='results', help='Results directory path')
    args = parser.parse_args()
    
    results_dir = args.dir
    
    print("=" * 80)
    print("AUTODOCK VINA VISUALIZATION")
    print("=" * 80)
    print(f"\nResults directory: {results_dir}")
    
    # Load data
    stats, combined = load_data(results_dir)
    
    if stats is None:
        print("ERROR: Could not load data.")
        return
    
    # Generate visualizations
    plot_boxplots(combined, results_dir)
    plot_heatmap(stats, results_dir)
    plot_bar_chart(stats, results_dir)
    save_detailed_table(stats, results_dir)
    
    print("\n" + "=" * 80)
    print("VISUALIZATION COMPLETE")
    print("=" * 80)
    print(f"\nGenerated files in {results_dir}:")
    print("  - plot1_binding_comparison.png")
    print("  - binding_affinity_matrix.png")
    print("  - plot3_validation_bars.png")
    print("  - detailed_comparison_table.csv")
    print("=" * 80)


if __name__ == "__main__":
    main()
