from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

base = Path(r'D:/PDAC_P1')
out = base / 'analysis/09_spatial_validation/SecAct'
summary = pd.read_csv(out / 'SecAct_five_targets_patient_summary_1000.tsv', sep='\t')
ecology = pd.read_csv(base / 'analysis/09_spatial_validation/GSE282302/spatial_ecology/GSE282302_spatial_ecology_patient.tsv', sep='\t')
targets = ['SPARCL1', 'VWF', 'MMRN2', 'ANGPT2', 'IL33']
summary = summary[summary['secreted'].isin(targets)].copy()
heat = summary.pivot(index='secreted', columns='dataset', values='positive_fraction').reindex(targets)
cellchat = pd.DataFrame({'candidate': targets, 'cellchat_edge': [0, 0, 0, 1, 0]})
cellchat = cellchat.set_index('candidate').reindex(targets)

plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams.update({'font.family': 'DejaVu Sans', 'axes.spines.top': False, 'axes.spines.right': False, 'figure.dpi': 180, 'savefig.dpi': 300})
fig = plt.figure(figsize=(10.5, 7.0))
grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1], hspace=0.42, wspace=0.28)

ax1 = fig.add_subplot(grid[0, 0])
im1 = ax1.imshow(heat.values, cmap='Blues', vmin=0, vmax=1, aspect='auto')
ax1.set_xticks(range(heat.shape[1]), heat.columns)
ax1.set_yticks(range(heat.shape[0]), heat.index)
for i in range(heat.shape[0]):
    for j in range(heat.shape[1]):
        ax1.text(j, i, f'{heat.iloc[i, j]:.2f}', ha='center', va='center', color='white' if heat.iloc[i, j] > 0.55 else '#222222')
fig.colorbar(im1, ax=ax1, fraction=0.046, pad=0.04, label='Patient/sample positive fraction')
ax1.set_xticks(range(heat.shape[1])); ax1.set_xticklabels(heat.columns)
ax1.set_yticks(range(heat.shape[0])); ax1.set_yticklabels(heat.index)
ax1.set_title('A  SecAct direction across spatial cohorts', loc='left', fontweight='bold')
ax1.set_xlabel('Cohort')
ax1.set_ylabel('Secreted protein')

ax2 = fig.add_subplot(grid[0, 1])
evidence = pd.DataFrame({'SecAct cross-cohort': [1] * len(targets), 'CellChat immune edge': [0, 0, 0, 1, 0], 'Immune neighborhood': [1] * len(targets)}, index=targets)
cmap = LinearSegmentedColormap.from_list('evidence', ['#E8EDF2', '#D55E00'])
im2 = ax2.imshow(evidence.values, cmap=cmap, vmin=0, vmax=1, aspect='auto')
ax2.set_xticks(range(evidence.shape[1]), evidence.columns)
ax2.set_yticks(range(evidence.shape[0]), evidence.index)
for i in range(evidence.shape[0]):
    for j in range(evidence.shape[1]):
        ax2.text(j, i, f'{evidence.iloc[i, j]:.0f}', ha='center', va='center')
ax2.set_xticks(range(evidence.shape[1])); ax2.set_xticklabels(['Cross-cohort\nSecAct', 'CellChat\nimmune edge', 'Immune\nneighborhood'])
ax2.set_yticks(range(evidence.shape[0])); ax2.set_yticklabels(evidence.index)
ax2.set_title('B  Integrated candidate prioritization', loc='left', fontweight='bold')
ax2.set_xlabel('Evidence layer')
ax2.tick_params(axis='x', labelrotation=0)
ax2.set_ylabel('')

ax3 = fig.add_subplot(grid[1, :])
plot = ecology[['median_myeloid_delta', 'median_T_NK_delta']].median().rename({'median_myeloid_delta': 'Myeloid', 'median_T_NK_delta': 'T/NK'}).reset_index()
plot.columns = ['program', 'delta']
bars = ax3.bar(plot['program'], plot['delta'], color=['#0072B2', '#009E73'], width=0.55)
ax3.axhline(0, color='#555555', linewidth=0.8)
ax3.set_ylabel('EPAS1-high minus EPAS1-low neighborhood score')
ax3.set_title('C  Immune neighborhood context in GSE282302', loc='left', fontweight='bold')
for bar, value in zip(bars, plot['delta']):
    ax3.text(bar.get_x() + bar.get_width() / 2, value + 0.003, f'{value:.3f}', ha='center', va='bottom', fontsize=10)
ax3.text(0.02, 0.92, '14/14 patients positive for both programs', transform=ax3.transAxes, fontsize=9)

fig.suptitle('EPAS1-associated spatial secreted signaling and immune ecology in PDAC', fontsize=14, fontweight='bold', y=0.98)
fig.savefig(out / 'Figure_integrated_EPAS1_SecAct_CellChat_ecology.pdf', bbox_inches='tight', pad_inches=0.15)
fig.savefig(out / 'Figure_integrated_EPAS1_SecAct_CellChat_ecology.png', bbox_inches='tight', pad_inches=0.15)
