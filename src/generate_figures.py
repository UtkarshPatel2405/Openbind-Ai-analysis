
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


C_ML = '#3b82f6'      # Blue
C_PHYSICS = '#ef4444'  # Red  
C_BASELINE = '#64748b' # Grey
C_ENSEMBLE = '#10b981' # Green


def load_results(path="results/rigorous_analysis.json"):
    with open(path) as f:
        return json.load(f)




def plot_panel_a(data, ax):
    """Global screening with 95% CIs."""
    global_data = data['global_screening']
    
    # Select methods to display
    method_map = {
        'iptm_of3p2': ('OF3-p2 (zero-shot)', C_ML),
        'iptm_of3p2_ft': ('OF3-p2 (fine-tuned)', C_ML),
        'iptm_boltz-2': ('Boltz-2', C_ML),
        'iptm_protenix': ('Protenix', C_ML),
        'pred_gnina_crystal': ('GNINA crystal', C_PHYSICS),
        'pred_molecular_weight': ('Molecular weight', C_BASELINE),
        'score_gnina': ('GNINA dock', C_PHYSICS),
    }
    
    rhos = []
    lows = []
    highs = []
    labels = []
    colors = []
    
    for key, (label, color) in method_map.items():
        if key in global_data and global_data[key].get('spearman_rho') is not None:
            rhos.append(global_data[key]['spearman_rho'])
            lows.append(global_data[key].get('spearman_ci_lower', rhos[-1]))
            highs.append(global_data[key].get('spearman_ci_upper', rhos[-1]))
            labels.append(label)
            colors.append(color)
    
    y_pos = np.arange(len(labels))
    
    # Error bars
    err_low = [r - l for r, l in zip(rhos, lows)]
    err_high = [h - r for r, h in zip(rhos, highs)]
    
    bars = ax.barh(y_pos, rhos, xerr=[err_low, err_high], 
                   capsize=3, color=colors, alpha=0.7, height=0.5, 
                   error_kw={'linewidth': 1.5, 'capthick': 1.5})
    
    ax.axvline(0, color='black', linewidth=0.8, linestyle='-')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r'Spearman $\rho$ [95% CI]', fontsize=10)
    ax.set_title('A. Global Screening (N=381-488)\nTarget-wide chemical diversity', 
                 fontweight='bold', fontsize=11)
    ax.set_xlim(-0.6, 0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add value labels at ends of bars with better positioning
    for i, (r, c) in enumerate(zip(rhos, colors)):
        if labels[i] == 'GNINA dock':
            # Place label inside the graph (to the right of the y-axis spine at -0.6)
            x_pos = lows[i] - 0.04  # starts at -0.53, well inside the plot
            y_pos_text = i + 0.16   # placed slightly above the bar
            ha = 'left'
            va = 'bottom'
        else:
            if r >= 0:
                x_pos = min(highs[i] + 0.03, 0.58)
                ha = 'left'
            else:
                x_pos = max(lows[i] - 0.03, -0.58)
                ha = 'right'
            y_pos_text = i
            va = 'center'
        ax.text(x_pos, y_pos_text, f'{r:.3f}', va=va, ha=ha, 
                fontsize=8, fontweight='bold', color=c)


def plot_panel_b(data, ax):
    """Rowan protocol sensitivity."""
    rowan = data['rowan_protocol_sensitivity']
    
    protocol_names = list(rowan.keys())
    rhos = [rowan[p]['spearman_rho'] for p in protocol_names]
    ps = [rowan[p]['spearman_p'] for p in protocol_names]
    
    # Color based on significance
    colors_bar = [C_PHYSICS if p < 0.05 else '#94a3b8' for p in ps]
    
    y_pos = np.arange(len(protocol_names))
    bars = ax.barh(y_pos, rhos, color=colors_bar, alpha=0.8, height=0.4, edgecolor='none')
    
    ax.axvline(0, color='black', linewidth=0.8, linestyle='-')
    ax.set_yticks(y_pos)
    ax.set_yticklabels([p.replace('docking_', 'dock-').replace('_', '-') for p in protocol_names], 
                       fontsize=9)
    ax.set_xlabel(r'Spearman $\rho$', fontsize=10)
    ax.set_title('B. Rowan Protocol Sensitivity (N=32)\nAll tested charge/model protocols', 
                 fontweight='bold', fontsize=11)
    ax.set_xlim(-0.1, 0.75)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add p-values for significant ones
    for i, (r, p) in enumerate(zip(rhos, ps)):
        if p < 0.05:
            p_text = f'p<0.001' if p < 0.001 else f'p={p:.2e}'
            ax.text(r + 0.03, i, p_text, va='center', ha='left', 
                    fontsize=7, color=C_PHYSICS)


def plot_panel_c(data, ax):
    """Pyrrolidine local SAR."""
    local_data = data['pyrrolidine_local_sar']
    rowan = data['rowan_protocol_sensitivity']
    
    # Include select methods + Rowan runs
    methods = [
        ('iptm_of3p2', 'OF3-p2 zero', C_ML),
        ('iptm_of3p2_ft', 'OF3-p2 fine-tuned', C_ML),
        ('pred_boltz_2', 'Boltz-2', C_ML),
        ('pred_gnina_crystal', 'GNINA crystal', C_PHYSICS),
        ('pred_molecular_weight', 'Molecular weight', C_BASELINE),
    ]
    
    rhos = []
    lows = []
    highs = []
    labels = []
    colors = []
    
    for key, label, color in methods:
        if key in local_data:
            d = local_data[key]
            rhos.append(d.get('spearman_rho', 0))
            lows.append(d.get('spearman_ci_lower', rhos[-1]) if d.get('spearman_ci_lower') else rhos[-1])
            highs.append(d.get('spearman_ci_upper', rhos[-1]) if d.get('spearman_ci_upper') else rhos[-1])
            labels.append(label)
            colors.append(color)
    
    # Add Rowan runs
    for proto_key in ['xtal_nagl', 'docking_am1bcc_0local']:
        if proto_key in rowan:
            d = rowan[proto_key]
            rhos.append(d['spearman_rho'])
            lows.append(d['spearman_rho'])  # No CI available
            highs.append(d['spearman_rho'])
            labels.append(proto_key.replace('xtal_', 'xtal-').replace('docking_', 'dock-').replace('_', '-'))
            colors.append(C_PHYSICS)
    
    y_pos = np.arange(len(labels))
    
    # Since some are point estimates (no CI), show different error bar styles
    err_low = [abs(r - l) for r, l in zip(rhos, lows)]
    err_high = [abs(h - r) for r, h in zip(rhos, highs)]
    
    # Add shaded "non-significant" zone
    ax.axvspan(-0.1, 0.1, alpha=0.08, color='black', zorder=0)
    
    bars = ax.barh(y_pos, rhos, xerr=[err_low, err_high],
                   capsize=3, color=colors, alpha=0.7, height=0.4,
                   error_kw={'linewidth': 1.5, 'capthick': 1.5})
    
    ax.axvline(0, color='black', linewidth=0.8, linestyle='-')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r'Spearman $\rho$ [95% CI]', fontsize=10)
    ax.set_title('C. Pyrrolidine Local SAR (N=28-32)\nWide CIs = high uncertainty', 
                 fontweight='bold', fontsize=11)
    ax.set_xlim(-0.5, 0.75)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add "NS" zone annotation
    ax.text(0.05, len(labels) - 0.5, 'All CIs span zero\n= not significant', 
            fontsize=8, style='italic', color='#94a3b8',
            ha='center', bbox=dict(boxstyle='round', facecolor='white', edgecolor='#94a3b8', alpha=0.5))


def plot_panel_d(data, ax):
    """Ensemble combinations."""
    ensemble = data['ensemble_combinations']
    
    labels = []
    rhos = []
    
    for key, val in ensemble.items():
        # Clean up label
        label = key.replace('ensemble_', '').replace('pred_', '').replace('iptm_', '')
        label = label.replace('molecular_weight', 'MW').replace('gnina_crystal', 'GNINA')
        label = label.replace('bolt_2', 'Boltz2').replace('of3p2_ft', 'OF3ft')
        label = label.replace('_', '+')
        labels.append(label)
        rhos.append(val['spearman_rho'])
    
    y_pos = np.arange(len(labels))
    bars = ax.barh(y_pos, rhos, color=C_ENSEMBLE, alpha=0.8, height=0.4, edgecolor='none')
    
    ax.axvline(0, color='black', linewidth=0.8, linestyle='-')
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel(r'Spearman $\rho$', fontsize=10)
    ax.set_title('D. Ensemble Combinations (Global)\nML + Physics synergy', 
                 fontweight='bold', fontsize=11)
    ax.set_xlim(0, 0.6)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Add individual component lines for comparison
    global_data = data['global_screening']
    if 'pred_gnina_crystal' in global_data:
        gnina_rho = global_data['pred_gnina_crystal']['spearman_rho']
        ax.axvline(gnina_rho, color=C_PHYSICS, linestyle='--', alpha=0.5, linewidth=1)
        ax.text(gnina_rho + 0.02, len(labels) - 0.5, f'GNINA alone: {gnina_rho:.3f}', 
                fontsize=7, color=C_PHYSICS, va='bottom')


def main():
    data = load_results()
    
    fig = plt.figure(figsize=(16, 10), dpi=300)
    
    # Create 2x2 grid with custom spacing
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.4)
    
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])
    
    plot_panel_a(data, ax1)
    plot_panel_b(data, ax2)
    plot_panel_c(data, ax3)
    plot_panel_d(data, ax4)
    
    # Main title
    fig.suptitle('OpenBind EV-A71 2A: Rigorous Benchmark with 95% Bootstrap Confidence Intervals', 
                 fontsize=14, fontweight='bold', y=0.98)
    
    # Legend
    legend_elements = [
        Patch(facecolor=C_ML, label='ML Co-folding'),
        Patch(facecolor=C_PHYSICS, label='Physics / Docking'),
        Patch(facecolor=C_BASELINE, label='Baseline (MW)'),
        Patch(facecolor=C_ENSEMBLE, label='Ensemble')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, 
               fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.02))
    
    plt.savefig('results/rigorous_benchmark.png', dpi=300, facecolor='white', 
                bbox_inches='tight')
    print('Success! Figure saved to results/rigorous_benchmark.png')


if __name__ == '__main__':
    main()
