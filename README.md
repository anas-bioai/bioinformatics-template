# Mitochondrial & Oxidative Stress Gene Dysregulation in Breast Cancer

## Overview
This project performs a complete RNA-seq differential expression analysis 
comparing breast cancer tumor tissue versus normal tissue, with a specific 
focus on mitochondrial and oxidative stress (ROS) gene dysregulation.

Developed as part of the Master in Computational Biology & Bioinformatics  
oral examination — FMPR.

---

## Biological Question
Are mitochondrial and oxidative stress genes differentially expressed 
in breast cancer compared to normal tissue?  
Is the observed pattern consistent with the **Warburg effect**?

---

## Dataset
| Parameter | Value |
|---|---|
| Normal samples | 14 |
| Tumor samples | 11 |
| Total genes | 58,037 |
| Source | RNA-seq count data (CSV format) |

---

## Pipeline
Raw counts → Log2 normalization → T-test → FDR correction →
Fold change → Gene annotation → Pathway filtering → Visualization
### Steps:
1. **Data Loading** — RNA-seq count files for Normal and Tumor samples
2. **Normalization** — log2(x+1) transformation
3. **Differential Expression** — independent t-test per gene
4. **FDR Correction** — Benjamini-Hochberg method
5. **Gene Annotation** — Ensembl ID → Gene Symbol via MyGene.info
6. **Pathway Filtering** — Mitochondrial & ROS gene families
7. **Visualization** — Volcano plots + Heatmap

---

## Key Results
| Finding | Value |
|---|---|
| Significant DEGs (FDR<0.05, \|log2FC\|>1) | 1,288 genes |
| Mitochondrial genes dysregulated | Downregulated in tumor |
| Key downregulated genes | ATP5PFP2, TIMM17B, COX7B2, MT-ND6 |
| Key ROS genes | CAT↑, NOX3↑, PRDX3P1↓ |
| Biological conclusion | Consistent with Warburg effect |

---

## Project Structure
```
mito-ros-analysis/
├── data/
│   ├── Normal/          ← normal tissue RNA-seq CSV files
│   └── Tumor/           ← tumor tissue RNA-seq CSV files
├── results/
│   ├── figures/
│   │   ├── volcano_overview.png     ← Volcano Plot Overview
│   │   ├── volcano_pathway.png      ← Volcano Plot Pathway Focus
│   │   └── heatmap.png              ← Heatmap Top 30 genes
│   ├── DEG_results.csv          ← full differential expression results
│   ├── significant_genes.csv    ← significant DEGs
│   ├── mito_genes.csv           ← mitochondrial genes
│   ├── ros_genes.csv            ← oxidative stress genes
│   └── results_annotated.csv    ← annotated results (offline cache)
└── analysis.ipynb       ← main analysis notebook
```
---

## Requirements

pandas
numpy
scipy
statsmodels
matplotlib
seaborn
mygene
adjustText

Install all dependencies:
```bash
pip install pandas numpy scipy statsmodels matplotlib seaborn mygene adjustText
```

---

## How to Run
1. Clone the repository
2. Place RNA-seq CSV files in `data/Normal/` and `data/Tumor/`
3. Open `analysis.ipynb` in Jupyter
4. Run all cells sequentially

> **Note:** Gene annotation (Section 5) requires internet connection  
> on first run. Results are cached in `results/results_annotated.csv`  
> for offline use after that.

---

## Biological Conclusion
Breast cancer tumor samples show systematic suppression of mitochondrial 
oxidative phosphorylation genes (ETC complexes I, IV, V) while some 
oxidative stress response genes remain active or are upregulated.  
This pattern is consistent with the **Warburg effect** — the metabolic 
reprogramming observed in cancer cells where mitochondrial respiration 
is reduced in favor of aerobic glycolysis.

---

## Author
**Anas**  
Master Computational Biology & Bioinformatics  
FMPR — 2025/2026
