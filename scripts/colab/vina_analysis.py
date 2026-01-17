"""
AutoDock Vina Results Analysis
==============================
Aggregates docking results and generates summary statistics.
Run after docking is complete.

Usage: python vina_analysis.py --dir /path/to/results
"""

import pandas as pd
import glob
import os
import argparse
from datetime import datetime


def format_number(num):
    """Format large numbers with commas"""
    return f"{num:,}"


def extract_affinities_from_logs(results_dir):
    """Extract binding affinities from Vina log files"""
    print("Extracting binding affinities from log files...")
    
    aggregated_dir = os.path.join(results_dir, 'aggregated')
    os.makedirs(aggregated_dir, exist_ok=True)
    
    raw_outputs = os.path.join(results_dir, 'raw_outputs')
    
    if not os.path.exists(raw_outputs):
        print(f"ERROR: raw_outputs directory not found at {raw_outputs}")
        return False
    
    # Process each ligand-pocket combination
    combinations = [d for d in os.listdir(raw_outputs) if os.path.isdir(os.path.join(raw_outputs, d))]
    
    for combo in combinations:
        output_dir = os.path.join(raw_outputs, combo)
        summary_file = os.path.join(aggregated_dir, f"{combo}_summary.csv")
        
        log_files = glob.glob(os.path.join(output_dir, "run_*.log"))
        if not log_files:
            continue
        
        with open(summary_file, 'w') as out:
            out.write("Run_Number,Mode,Affinity_kcal_mol,RMSD_lb,RMSD_ub\n")
            
            for log_file in log_files:
                run_num = os.path.basename(log_file).replace('run_', '').replace('.log', '')
                
                try:
                    with open(log_file) as f:
                        for line in f:
                            line = line.strip()
                            if line and line[0].isdigit() and '   ' in line:
                                parts = line.split()
                                if len(parts) >= 4:
                                    out.write(f"{run_num},{parts[0]},{parts[1]},{parts[2]},{parts[3]}\n")
                except Exception as e:
                    print(f"  Warning: Error reading {log_file}: {e}")
    
    print(f"✓ Extracted binding data to {aggregated_dir}/")
    return True


def load_and_analyze(results_dir):
    """Load summary files and generate statistics"""
    aggregated_dir = os.path.join(results_dir, 'aggregated')
    summary_files = glob.glob(os.path.join(aggregated_dir, '*_summary.csv'))
    
    if not summary_files:
        print("⚠️  No summary files found. Run extraction first.")
        return None
    
    all_results = []
    
    for f in summary_files:
        df = pd.read_csv(f)
        basename = os.path.basename(f).replace('_summary.csv', '')
        
        # Parse ligand and pocket from filename
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
    print(f"✓ Loaded {format_number(len(combined))} binding poses from {len(summary_files)} files")
    
    # Calculate statistics
    summary = combined.groupby(['Ligand', 'Pocket', 'Pocket_Name'])['Affinity_kcal_mol'].agg([
        ('Best_Affinity', 'min'),
        ('Mean_Affinity', 'mean'),
        ('Median_Affinity', 'median'),
        ('Std_Affinity', 'std'),
        ('Total_Poses', 'count')
    ]).reset_index()
    
    summary = summary.sort_values(['Ligand', 'Best_Affinity'])
    
    return summary, combined


def validate_results(summary):
    """Check if P0 (Orthosteric) shows best binding as expected"""
    print("\n" + "=" * 80)
    print("                     VALIDATION CHECK")
    print("=" * 80)
    
    validation_passed = True
    
    for ligand in summary['Ligand'].unique():
        ligand_data = summary[summary['Ligand'] == ligand]
        best_pocket = ligand_data.loc[ligand_data['Best_Affinity'].idxmin(), 'Pocket']
        best_affinity = ligand_data['Best_Affinity'].min()
        
        p0_data = ligand_data[ligand_data['Pocket'] == 'P0']
        p0_affinity = p0_data['Best_Affinity'].values[0] if len(p0_data) > 0 else None
        
        print(f"\n{ligand}:")
        print(f"  Best binding at: {best_pocket}")
        print(f"  Best affinity: {best_affinity:.2f} kcal/mol")
        
        if best_pocket == 'P0':
            print(f"  ✓ VALIDATION PASSED: P0 (Orthosteric) shows best binding")
        else:
            validation_passed = False
            if p0_affinity:
                delta = p0_affinity - best_affinity
                print(f"  ⚠️  WARNING: {best_pocket} shows better binding than P0")
                print(f"     {best_pocket}: {best_affinity:.2f} kcal/mol")
                print(f"     P0: {p0_affinity:.2f} kcal/mol")
                print(f"     Difference: {delta:.2f} kcal/mol")
    
    return validation_passed


def save_summary(summary, results_dir):
    """Save final summary to CSV"""
    output_file = os.path.join(results_dir, 'FINAL_SUMMARY.csv')
    summary.to_csv(output_file, index=False)
    print(f"\n✓ Summary saved to: {output_file}")
    return output_file


def generate_report(summary, results_dir):
    """Generate detailed text report"""
    report_file = os.path.join(results_dir, 'analysis', 'statistical_report.txt')
    os.makedirs(os.path.dirname(report_file), exist_ok=True)
    
    with open(report_file, 'w') as f:
        f.write("=" * 80 + "\n")
        f.write("AutoDock Vina - Phase 1 Validation Report\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        for ligand in summary['Ligand'].unique():
            ligand_data = summary[summary['Ligand'] == ligand].sort_values('Best_Affinity')
            
            f.write(f"\n{'='*80}\n")
            f.write(f"LIGAND: {ligand}\n")
            f.write(f"{'='*80}\n\n")
            
            f.write("Ranking by Best Affinity:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Rank':<6} {'Pocket':<8} {'Site':<25} {'Best':<10} {'Mean':<10} {'Std':<8}\n")
            f.write("-" * 80 + "\n")
            
            for rank, (idx, row) in enumerate(ligand_data.iterrows(), 1):
                f.write(f"{rank:<6} {row['Pocket']:<8} {row['Pocket_Name']:<25} "
                       f"{row['Best_Affinity']:>8.2f}  {row['Mean_Affinity']:>8.2f}  "
                       f"{row['Std_Affinity']:>6.2f}\n")
            
            f.write("\n")
    
    print(f"✓ Report saved to: {report_file}")
    return report_file


def main():
    parser = argparse.ArgumentParser(description='Analyze AutoDock Vina results')
    parser.add_argument('--dir', default='results', help='Results directory path')
    args = parser.parse_args()
    
    results_dir = args.dir
    
    print("=" * 80)
    print("AUTODOCK VINA RESULTS ANALYSIS")
    print("=" * 80)
    print(f"\nResults directory: {results_dir}")
    
    # Step 1: Extract affinities
    extract_affinities_from_logs(results_dir)
    
    # Step 2: Load and analyze
    result = load_and_analyze(results_dir)
    if result is None:
        return
    
    summary, combined = result
    
    # Step 3: Print summary
    print("\n" + "=" * 80)
    print("                     VALIDATION SUMMARY")
    print("=" * 80)
    print(summary.to_string(index=False))
    
    # Step 4: Validate
    validate_results(summary)
    
    # Step 5: Save outputs
    save_summary(summary, results_dir)
    generate_report(summary, results_dir)
    
    print("\n" + "=" * 80)
    print("✓ Analysis complete!")
    print("=" * 80)


if __name__ == "__main__":
    main()
