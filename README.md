# Bulk-RNAseq-Humanized-ImmunoRT

This repository contains scripts and functions used to process and analyze the data featured in the publication:

📄 **Cogels M. et al., [Publication Title] ([Link])**

---

## 📁 Repository Structure

### R
Contains R scripts for:
- **Preprocessing & Normalization**: Cleaning and formatting count and TPM files from human sequencing data (`clean_human_count_TPM.R`).
- **Differential Expression & Enrichment Analysis**: Performing Differential Gene Expression Analysis (DGEA) and Gene Set Enrichment Analysis (GSEA) (`DGEA_and_GSEA.R`).
- **GSVA Analysis**: Conducting Gene Set Variation Analysis (GSVA) (`GSVA.R`).
- **xCell & Heatmaps**: Performing xCell analysis and generating heatmaps to visualize xCell results and specific gene sets (`xcell_genesets_heatmaps.R`).

### input_files
Contains:
- **Gene Sets File**: Excel file with gene sets used in `xcell_genesets_heatmaps.R` (`25_02_05_RNAseq_Final_Heatmaps_Clean_V2.xlsx`).

---

## 📊 Data Availability
All additional data required to replicate the analyses can be accessed at:
🔗 **[Zenodo / Data Repository Link]**

---

## 🛠️ Software Requirements
The scripts in this repository have been tested with the following software versions:
- **R**: v4.4.0

---

## 📢 Citation
If you use these scripts and analyses in your research, please cite:

📄 **Cogels M., et al., [Title of the Paper]**
📌 **Publication Link: [Link]**

---

## 📘 Further Information
This analysis was conducted by **Matteo Serra**, **David Venet** and **Morgane Cogels** as part of research at **Institut Jules Bordet**.

For inquiries, feel free to contact 📩 **matteo95serra@gmail.com** or **morgane.cogels@gmail.com** or open an issue on this repository.

