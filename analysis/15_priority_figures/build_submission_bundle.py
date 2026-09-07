from pathlib import Path
import shutil
import matplotlib.pyplot as plt
from matplotlib.image import imread

ROOT = Path(r"D:/第二篇大论文")
OUT = ROOT / "analysis/15_priority_figures/submission_bundle_2026-09-06_v3"
SUPP = OUT / "supplementary_figures"
MAIN = OUT / "main_figures"
SUPP.mkdir(parents=True, exist_ok=True)
MAIN.mkdir(parents=True, exist_ok=True)

main_src = {
    1: ROOT / "analysis/15_priority_figures/main_figures_relayout/Figure_1_scRNA_discovery_updated.png",
    2: ROOT / "analysis/15_priority_figures/main_figures_unified/Figure_2_bulk_spatial_validation_unified.png",
    3: ROOT / "analysis/15_priority_figures/main_figures_unified/Figure_3_secact_immune_ecology_unified.png",
    4: ROOT / "analysis/21_revision_strengthening/Figure_4_virtual_perturbation_unified_REVISED.pdf",
    5: ROOT / "analysis/21_revision_strengthening/Figure_5_drug_structure_closure_unified_REVISED.pdf",
    6: ROOT / "analysis/15_priority_figures/main_figures_unified/Figure_6_EPAS1_closed_loop_unified.png",
}

for number, path in main_src.items():
    if not path.exists():
        raise FileNotFoundError(path)
    shutil.copy2(path, MAIN / f"Figure_{number}{path.suffix}")

src = {
    1: ROOT / "analysis/15_priority_figures/Figure_S6_scRNA_discovery_overview.pdf",
    2: ROOT / "analysis/09_spatial_validation/SecAct/Figure_integrated_EPAS1_SecAct_CellChat_ecology.pdf",
    3: ROOT / "analysis/08_program_validation/GSE21501/GSE21501_EPAS1_external_survival_KM.pdf",
    4: ROOT / "analysis/15_priority_figures/Figure_S11_CPTAC_protein_evidence.pdf",
    5: ROOT / "analysis/15_priority_figures/Figure_S13_negative_pharmacology_gates.pdf",
    6: ROOT / "analysis/15_priority_figures/Figure_S14_survival_ML_sensitivity.pdf",
    7: ROOT / "analysis/15_priority_figures/Figure_S16_EPAS1_structure_docking_QC.pdf",
    8: ROOT / "analysis/15_priority_figures/Figure_S17_MR_instrument_gate.pdf",
    9: ROOT / "analysis/15_priority_figures/Figure_S18_GSE21501_survival_forest.pdf",
    10: ROOT / "analysis/15_priority_figures/Figure_S19_functional_immune_WGCNA.pdf",
    11: ROOT / "analysis/15_priority_figures/Figure_S20_TCGA_bulk_QC_PCA.pdf",
    12: ROOT / "analysis/15_priority_figures/Figure_S21_spatial_image_montage.pdf",
    13: ROOT / "analysis/15_priority_figures/Figure_S22_endothelial_DE_pathway_gate.pdf",
    14: ROOT / "analysis/15_priority_figures/Figure_S23_EPAS1_complex_binding_modes.pdf",
    15: ROOT / "analysis/15_priority_figures/Figure_S15_scRNA_annotation_CopyKAT_CNV.pdf",
    16: ROOT / "analysis/15_priority_figures/scrna_additional/Figure_S15_scRNA_sample_QC.pdf",
    17: ROOT / "analysis/15_priority_figures/scrna_additional/Figure_S16_scRNA_marker_DotPlot.pdf",
    18: ROOT / "analysis/15_priority_figures/scrna_additional/Figure_S17_endothelial_pseudobulk_volcano.pdf",
    19: ROOT / "analysis/15_priority_figures/Figure_S7_endothelial_candidate_convergence.pdf",
    20: ROOT / "analysis/15_priority_figures/Figure_S8_bulk_program_validation.pdf",
    21: ROOT / "analysis/15_priority_figures/Figure_S9_spatial_ecology_patient_level.pdf",
    22: ROOT / "analysis/15_priority_figures/Figure_S10_secact_angpt2_replication.pdf",
    23: ROOT / "analysis/15_priority_figures/Figure_S12_DrugReflector_prioritization.pdf",
    24: ROOT / "analysis/15_priority_figures/Figure_S16_EPAS1_main_docking_QC.pdf",
    25: ROOT / "analysis/15_priority_figures/figure2_supplement/Figure_S24_bulk_QC_DE_GSEA.pdf",
    26: ROOT / "analysis/15_priority_figures/figure2_supplement/Figure_S25_spatial_domains_colocalization.pdf",
    28: ROOT / "analysis/18_external_strengthening/External_Figure_GSE202051_endothelial_localization.pdf",
    29: ROOT / "analysis/19_translation_closure/Supplementary_Figure_ANGPT2_cross_method_evidence.pdf",
    30: ROOT / "analysis/20_GSE205049_immune_receiver/Supplementary_Figure_GSE205049_immune_receiver_sensitivity.pdf",
    31: ROOT / "analysis/21_revision_strengthening/Figure_S31_EPAS1_directionality_audit.pdf",
    32: ROOT / "analysis/21_revision_strengthening/Figure_S32_spatial_conditional_model.pdf",
}

for number, path in src.items():
    if not path.exists():
        raise FileNotFoundError(path)
    shutil.copy2(path, SUPP / f"Supplementary_Figure_S{number:02d}.pdf")
    for ext in (".png", ".tiff", ".tif"):
        candidate = path.with_suffix(ext)
        if candidate.exists():
            shutil.copy2(candidate, SUPP / f"Supplementary_Figure_S{number:02d}{ext}")

# S27 combines the two complementary virtual-perturbation audits into one
# upload item while retaining the original source files outside the bundle.
ko_png = ROOT / "analysis/13_virtual_perturbation/closed_loop/supplement_gene_level/Figure_S26_EPAS1_KO_gene_level_audit.png"
net_png = ROOT / "analysis/13_virtual_perturbation/closed_loop/supplement_gene_level/Figure_S27_EPAS1_plotKO_networks.png"
fig = plt.figure(figsize=(12, 15), facecolor="white")
fig.text(0.5, 0.985, "Supplementary Figure S27. EPAS1 virtual-perturbation gene and network audit", ha="center", va="top", fontsize=15, fontweight="bold")
ax1 = fig.add_axes([0.04, 0.42, 0.92, 0.53]); ax1.imshow(imread(ko_png)); ax1.axis("off"); ax1.text(0.0, 1.01, "A", transform=ax1.transAxes, fontsize=14, fontweight="bold", va="bottom")
ax2 = fig.add_axes([0.04, 0.04, 0.92, 0.30]); ax2.imshow(imread(net_png)); ax2.axis("off"); ax2.text(0.0, 1.01, "B", transform=ax2.transAxes, fontsize=14, fontweight="bold", va="bottom")
fig.savefig(SUPP / "Supplementary_Figure_S27.png", dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(SUPP / "Supplementary_Figure_S27.pdf", bbox_inches="tight", facecolor="white")
fig.savefig(SUPP / "Supplementary_Figure_S27.tiff", dpi=600, bbox_inches="tight", facecolor="white")
plt.close(fig)

(OUT / "README_upload_manifest.txt").write_text(
    "Canonical upload sequence: Supplementary_Figure_S01 through Supplementary_Figure_S32.\n"
    "Main figures are in main_figures/Figure_1 through Figure_6; revised Figures 4 and 5 are included.\n"
    "PDF is the primary upload candidate; PNG/TIFF are high-resolution delivery copies.\n"
    "Original legacy files remain in their analysis directories and are not upload assets.\n",
    encoding="utf-8",
)
