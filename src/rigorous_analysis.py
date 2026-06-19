"""
Rigorous Benchmark Analysis for OpenBind EV-A71 2A Protease
============================================================

Performs:
1. Bootstrap confidence intervals for all correlations
2. Systematic protocol sensitivity (all Rowan runs)
3. Ensemble analysis: ML + physics combinations
4. Exports fully reproducible results

Run: python src/rigorous_analysis.py
"""
import os, json, warnings
import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import pearsonr, spearmanr

warnings.filterwarnings('ignore')
np.random.seed(42)

# Paths
OPENBIND_DIR = "EV-A71_2A_benchmark"
ROWAN_DIR = "openbind_ev_a71_rbfe"


def bootstrap_ci(x, y, stat_func, n_bootstrap=2000, ci=0.95):
    """Bootstrap confidence interval for correlation statistic."""
    n = len(x)
    if n < 10:
        return np.nan, np.nan, np.nan
    stats = []
    for _ in range(n_bootstrap):
        idx = np.random.choice(n, size=n, replace=True)
        xb, yb = x[idx], y[idx]
        if len(np.unique(xb)) > 1 and len(np.unique(yb)) > 1:
            try:
                s = stat_func(xb, yb)
                if not np.isnan(s):
                    stats.append(s)
            except:
                pass
    if len(stats) < 100:
        return np.nan, np.nan, np.nan
    stats_arr = np.array(stats)
    lo = np.percentile(stats_arr, (1 - ci) / 2 * 100)
    hi = np.percentile(stats_arr, (1 + ci) / 2 * 100)
    return stat_func(x, y), lo, hi


def corr_stats(x, y, n_boot=2000):
    """Compute Spearman + Pearson with p-values and bootstrap 95% CIs."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = ~(np.isnan(x) | np.isnan(y))
    x, y = x[mask], y[mask]
    N = len(x)
    
    if N < 10 or len(np.unique(x)) < 2 or len(np.unique(y)) < 2:
        return {"N": N}
    
    rho, p_s = spearmanr(x, y)
    r, p_p = pearsonr(x, y)
    
    _, rho_lo, rho_hi = bootstrap_ci(x, y, lambda a, b: spearmanr(a, b)[0])
    _, r_lo, r_hi = bootstrap_ci(x, y, lambda a, b: pearsonr(a, b)[0])
    
    return {
        "spearman_rho": float(rho), "spearman_p": float(p_s),
        "spearman_ci_lower": float(rho_lo) if not np.isnan(rho_lo) else None,
        "spearman_ci_upper": float(rho_hi) if not np.isnan(rho_hi) else None,
        "pearson_r": float(r), "pearson_p": float(p_p),
        "pearson_ci_lower": float(r_lo) if not np.isnan(r_lo) else None,
        "pearson_ci_upper": float(r_hi) if not np.isnan(r_hi) else None,
        "N": N
    }


def load_data():
    """Load and canonicalize all benchmark data."""
    print("[INFO] Loading benchmark data...")
    
    ref = pd.read_csv(os.path.join(OPENBIND_DIR, "affinity", "reference", "fragalysis_compound_reference.csv"))
    cofold = pd.read_parquet(os.path.join(OPENBIND_DIR, "structure", "processed_outputs", "final_cofolding_pose_data.parquet"))
    dock = pd.read_parquet(os.path.join(OPENBIND_DIR, "structure", "processed_outputs", "final_docking_pose_data.parquet"))
    preds = pd.read_csv(os.path.join(OPENBIND_DIR, "affinity", "outputs", "compound_level_prediction_analysis.csv"))
    
    pyr = pd.read_csv(os.path.join(ROWAN_DIR, "openbind_ev71_2a_pyrrolidine_benchmark_release", "subset", "pyrrolidine_32_subset.csv"))
    methyl = pd.read_csv(os.path.join(ROWAN_DIR, "openbind_ev71_2a_methylthio_pyrimidine_benchmark_release", "subset", "methylthio_pyrimidine_76_subset.csv"))
    
    smiles_cache = {}
    def canon(s):
        s_clean = s.split()[0]
        if s_clean not in smiles_cache:
            mol = Chem.MolFromSmiles(s_clean)
            smiles_cache[s_clean] = Chem.MolToSmiles(mol, canonical=True) if mol else s_clean
        return smiles_cache[s_clean]
    
    for df in [ref, cofold, dock, preds, pyr, methyl]:
        df['canon_smiles'] = (df['smiles'] if 'smiles' in df.columns else df['ligand_smiles']).apply(canon)
    
    return ref, cofold, dock, preds, pyr, methyl


def analyze_global(df_global):
    """Analyze global screening correlations."""
    print("\n[SECTION] Global Screening (N=381-488)")
    print("=" * 80)
    
    results = {}
    for col in df_global.columns:
        if col not in ['canon_smiles', 'experimental_pKD']:
            sub = df_global[['experimental_pKD', col]].dropna()
            if len(sub) > 10:
                res = corr_stats(sub['experimental_pKD'].values, sub[col].values)
                results[col] = res
                sig = "***" if res.get('spearman_p', 1) < 0.001 else "**" if res.get('spearman_p', 1) < 0.01 else "*" if res.get('spearman_p', 1) < 0.05 else ""
                ci = f"[{res.get('spearman_ci_lower', 'N/A'):.3f}, {res.get('spearman_ci_upper', 'N/A'):.3f}]"
                print(f"{col:<30} | N={res['N']:<3} | rho={res.get('spearman_rho', 'N/A'):>6.3f}{' ' + sig:<4} | CI={ci} | p={res.get('spearman_p', 'N/A'):.2e}")
    return results


def analyze_local(df_pyr):
    """Analyze pyrrolidine local SAR."""
    print("\n[SECTION] Pyrrolidine Local SAR (N=28-32)")
    print("=" * 80)
    
    results = {}
    for col in df_pyr.columns:
        if col not in ['canon_smiles', 'experimental_pKD']:
            sub = df_pyr[['experimental_pKD', col]].dropna()
            if len(sub) > 5:
                res = corr_stats(sub['experimental_pKD'].values, sub[col].values)
                results[col] = res
                sig = "***" if res.get('spearman_p', 1) < 0.001 else "**" if res.get('spearman_p', 1) < 0.01 else "*" if res.get('spearman_p', 1) < 0.05 else ""
                ci_l = res.get('spearman_ci_lower')
                ci_u = res.get('spearman_ci_upper')
                ci_str = f"[{ci_l:.3f}, {ci_u:.3f}]" if ci_l is not None and ci_u is not None else "[N/A, N/A]"
                print(f"{col:<30} | N={res['N']:<3} | rho={res.get('spearman_rho', 'N/A'):>6.3f}{' ' + sig:<4} | CI={ci_str} | p={res.get('spearman_p', 'N/A'):.2e}")
    return results


def analyze_rowan_protocols():
    """Load and report ALL Rowan run statistics."""
    print("\n[SECTION] Rowan Protocol Sensitivity (all runs, N=32)")
    print("=" * 80)
    
    rowan_csv = os.path.join(ROWAN_DIR, "openbind_ev71_2a_pyrrolidine_benchmark_release", "results", "rowan_results_overall_by_run.csv")
    df = pd.read_csv(rowan_csv)
    
    results = {}
    for _, r in df.iterrows():
        name = r['run']
        results[name] = {
            "spearman_rho": float(r['spearman_rho']), "spearman_p": float(r['spearman_p']),
            "pearson_r": float(r['pearson_r']), "pearson_p": float(r['pearson_p']),
            "r2_score": float(r['r2_score']), "N": int(r['n'])
        }
        sig = "***" if r['spearman_p'] < 0.001 else "**" if r['spearman_p'] < 0.01 else "*" if r['spearman_p'] < 0.05 else ""
        print(f"{name:<30} | N={int(r['n']):<3} | rho={r['spearman_rho']:>6.3f} {sig:<3} | R2={r['r2_score']:>6.3f} | p={r['spearman_p']:.2e}")
    return results


def analyze_ensemble(df_global, methods_to_test=None):
    """Test ensemble combinations of scoring methods."""
    print("\n[SECTION] Ensemble Analysis (Global)")
    print("=" * 80)
    
    # Candidate pairs: (ML_method, physics_method)
    pairs = methods_to_test or [
        ('pred_boltz_2', 'pred_gnina_crystal'),
        ('pred_molecular_weight', 'pred_boltz_2'),
        ('iptm_of3p2_ft', 'pred_gnina_crystal')
    ]
    
    results = {}
    for m1, m2 in pairs:
        if m1 in df_global.columns and m2 in df_global.columns:
            sub = df_global[['canon_smiles', 'experimental_pKD', m1, m2]].dropna()
            if len(sub) > 10:
                # Normalize both to z-scores
                z1 = (sub[m1] - sub[m1].mean()) / sub[m1].std()
                z2 = (sub[m2] - sub[m2].mean()) / sub[m2].std()
                ensemble = 0.5 * z1 + 0.5 * z2
                res = corr_stats(sub['experimental_pKD'].values, ensemble.values)
                label = f"ensemble_{m1}_{m2}"
                results[label] = res
                sig = "***" if res.get('spearman_p', 1) < 0.001 else "**" if res.get('spearman_p', 1) < 0.01 else "*" if res.get('spearman_p', 1) < 0.05 else ""
                ci_l = res.get('spearman_ci_lower')
                ci_u = res.get('spearman_ci_upper')
                ci_str = f"[{ci_l:.3f}, {ci_u:.3f}]" if ci_l is not None and ci_u is not None else "[N/A, N/A]"
                print(f"{label:<40} | N={res['N']:<3} | rho={res.get('spearman_rho', 'N/A'):>6.3f}{' ' + sig:<4} | CI={ci_str} | p={res.get('spearman_p', 'N/A'):.2e}")
    return results


def build_global_df(ref, cofold, dock, preds):
    """Build the global compound-level dataset."""
    df = ref.groupby('canon_smiles', as_index=False)['experimental_pKD'].mean()
    
    # Affinity predictions
    for m in preds['method'].unique():
        m_df = preds[preds['method'] == m].groupby('canon_smiles')['predicted_affinity'].mean().reset_index()
        m_df.columns = ['canon_smiles', f'pred_{m}']
        df = df.merge(m_df, on='canon_smiles', how='left')
    
    # Co-folding confidence scores (top-ranked pose only)
    cofold_r0 = cofold[cofold['rank'] == 0]
    for m in cofold_r0['method'].unique():
        m_df = cofold_r0[cofold_r0['method'] == m].groupby('canon_smiles')['rank_score'].mean().reset_index()
        m_df.columns = ['canon_smiles', f'iptm_{m}']
        df = df.merge(m_df, on='canon_smiles', how='left')
    
    # GNINA docking score
    dock_r0 = dock[dock['rank'] == 0]
    gnina = dock_r0[dock_r0['method'] == 'gnina_multi'].groupby('canon_smiles')['rank_score'].mean().reset_index()
    gnina.columns = ['canon_smiles', 'score_gnina']
    df = df.merge(gnina, on='canon_smiles', how='left')
    
    return df


def save_results(global_res, local_res, rowan_res, ensemble_res, output_dir="results"):
    """Save all results to JSON."""
    os.makedirs(output_dir, exist_ok=True)
    
    output = {
        "global_screening": global_res,
        "pyrrolidine_local_sar": local_res,
        "rowan_protocol_sensitivity": rowan_res,
        "ensemble_combinations": ensemble_res
    }
    
    with open(os.path.join(output_dir, "rigorous_analysis.json"), "w") as f:
        json.dump(output, f, indent=2)
    
    print(f"\n[INFO] Results saved to {output_dir}/rigorous_analysis.json")
    return output


def main():
    ref, cofold, dock, preds, pyr, methyl = load_data()
    
    # Build global dataset
    df_global = build_global_df(ref, cofold, dock, preds)
    
    # Run all analyses
    global_res = analyze_global(df_global)
    
    df_pyr = df_global.merge(pyr[['canon_smiles']], on='canon_smiles')
    local_res = analyze_local(df_pyr)
    
    rowan_res = analyze_rowan_protocols()
    
    ensemble_res = analyze_ensemble(df_global)
    
    # Save
    output = save_results(global_res, local_res, rowan_res, ensemble_res)
    
    print("\n[SUMMARY]")
    print(f"  Global methods analyzed: {len(global_res)}")
    print(f"  Local (pyrrolidine) methods: {len(local_res)}")
    print(f"  Rowan protocols: {len(rowan_res)}")
    print(f"  Ensemble combinations: {len(ensemble_res)}")
    print("\n[INFO] Next: run src/generate_figures.py to create publication-quality plots")


if __name__ == "__main__":
    main()
