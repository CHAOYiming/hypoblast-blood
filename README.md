# Extraembryonic mesoderm from the hypoblast facilitates human pregastrula development
Early human embryos grow before circulation is established, but whether extraembryonic tissues provide interim blood-like support is unknown. Here, using human stem-cell-based embryo models, we identify an extraembryonic mesoderm at pregastrula that expresses embryonic hemoglobin. Molecular barcoding, direct differentiation of hypoblast, and live-imaging identify that extraembryonic mesoderm emerges from the hypoblast lineage, which alleviates hypoxic stress. CDX2 marks this transition, and its loss reduces hemoglobin expression and increases hypoxic stress in embryo models. CRISPR activation screening identifies an LMO2 regulatory network that confers on hypoblasts the erythroid programme. Human embryo transcriptomics and immunostaining detected this population in vivo. These findings reveal a hypoblast-derived extraembryonic mesoderm that confers blood-like function in pregastrula before circulatory maturation.

## The download link of processed files: 
Human embryo scRNA-seq dataset [PMID: 39543283] can be downloaded from Fredrik Lanner’s lab website: 
https://petropoulos-lanner-labs.clintec.ki.se/dataset.download.html

Mouse embryo dataset [PMID: 30787436] can be retrieved from the MouseGastrulationData repository:
https://github.com/MarioniLab/MouseGastrulationData

Human blastocyst dataset [PMID: 38277271] can be downloaded from GEO: 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE136106

The in-house sequencing matrix has been deposited here: 
1. HCEB from d1-d10: https://www.dropbox.com/scl/fo/58bih6gikizyh8exh9agh/AGShN5qq10-mmSc4Ar-srHA?rlkey=n1xdqcwh7yckchh9me521lfud&st=y1p0ubi0&dl=0.
2. larry-barcoded HCEB from two timepoints (d4, d12): https://drive.google.com/drive/folders/1oQZTiNSTdNOoNjWJvlpcrbJNxnDnUBP-?usp=drive_link
3. WT and CDX2-KO Peri-gastruloid from two timepoints (WT: d4/d8; KO: d6/d8): https://drive.google.com/drive/folders/13x08RM6DQTaOyLVfwWwPyNl0xN3A5jKi?usp=drive_link

Backup: https://doi.org/10.6084/m9.figshare.32979542

## Analysis scripts in this repository

The repository contains analysis workflows for public embryo datasets, human
stem-cell-derived embryo models, LARRY lineage tracing, peri-gastruloids,
Waddington-OT analysis, and the CRISPR activation screen.

```text
hypoblast-blood/
│
├── public_data/
│   ├── 01_human_embryo_reanalysis.R
│   ├── 02_mouse_embryo_reanalysis.R
│   ├── 03_icm_bulk.R
│   ├── 04_lcxyw_construction.ipynb
│   ├── 05_lcxyw_hceb_exm_mapping.ipynb
│   ├── 06_cs6_hceb_exm_mapping.ipynb
│   └── 07_lanner_wot_full_visualization.ipynb
│
├── HCEB/
│   ├── 01_qc.ipynb
│   ├── 02_integration.ipynb
│   ├── 03_annotation.ipynb
│   ├── 04_blood_trajectory.ipynb
│   ├── 05_endo_subtype.ipynb
│   ├── 06_tvae_stage_d4_10.ipynb
│   └── 03_hceb_exm_visualization.ipynb
│
├── LARRY/
│   ├── 01_qc.ipynb
│   ├── 02_annotation.ipynb
│   ├── 03_barcode_preprocess.ipynb
│   ├── 04_clone_analysis.ipynb
│   ├── 05_lineage.py
│   ├── 01_hceb_larry_annotation_visualization.ipynb
│   ├── 03_hceb_larry_basic_clonality.ipynb
│   └── 04_hceb_larry_clone_visualization.ipynb
│
├── PeriGas/
│   ├── 01_qc.ipynb
│   ├── 02_integration.ipynb
│   ├── 03_annotation.ipynb
│   ├── 04_wt_ko.ipynb
│   ├── 05_function.ipynb
│   └── 06_tvae_stage_wt.ipynb
│
└── CRISPRa/
    ├── 01_CRISPRa_preprocess.sh
    ├── 02_CRISPRa_downstream_analysis.R
    └── 03_CRISPRa_screen_library_summary.xls
```

## Software and package versions

The main computational analyses were performed using:

- Python 3.10.20
- Scanpy 1.11.5
- AnnData 0.11.4
- NumPy 2.2.6
- pandas 2.3.3
- SciPy 1.15.3
- scikit-learn 1.7.2
- Matplotlib 3.10.9
- BBKNN 1.6.0
- leidenalg 0.11.0
- umap-learn 0.5.12
- Waddington-OT (`wot`) 1.0.8.post2
- POT 0.9.6.post1
- statsmodels 0.14.6

The full Python/conda environment is provided in `environment.yml`.

## Figure-to-script mapping

| Manuscript panel | Analysis | Script / notebook |
|---|---|---|
| Fig. 1E–F | HCEB single-cell integration, annotation, and characterization of extraembryonic mesoderm states | `HCEB/02_integration.ipynb`; `HCEB/03_annotation.ipynb`; `HCEB/03_hceb_exm_visualization.ipynb` |
| Extended Data Fig. 1A | Marker-gene visualization for HCEB cell-type annotation | `HCEB/03_annotation.ipynb`; `HCEB/03_hceb_exm_visualization.ipynb` |
| Extended Data Fig. 1B | Cell-type composition across the HCEB time course | `HCEB/03_hceb_exm_visualization.ipynb` |
| Extended Data Fig. 1C | TemporalVAE developmental-stage mapping of HCEBs | `HCEB/06_tvae_stage_d4_10.ipynb` |
| Fig. 2A–B | Construction and visualization of the integrated human embryo reference atlas | `public_data/04_lcxyw_construction.ipynb` |
| Fig. 2C–F | Projection of HCEB-derived ExM populations onto the integrated human embryo reference and ExM-program comparison | `public_data/05_lcxyw_hceb_exm_mapping.ipynb` |
| Fig. 2G–H | Reanalysis of the CS6 human embryo and projection of HCEB-derived Trans ExM | `public_data/06_cs6_hceb_exm_mapping.ipynb` |
| Fig. 2I | Slingshot trajectory inference across epiblast-, hypoblast-, and mesoderm-associated populations | `public_data/07_lanner_slingshot.R` |
| Fig. 2J | Waddington-OT fate-probability analysis of early hypoblast | `public_data/07_lanner_wot_full_visualization.ipynb` |
| Extended Data Fig. 2A | Early human embryo lineage-reference visualization | `public_data/04_lcxyw_construction.ipynb` |
| Extended Data Fig. 2B–C | Slingshot lineage trajectories and lineage summaries | `public_data/07_lanner_slingshot.R` |
| Extended Data Fig. 2D | Waddington-OT fate-probability analysis | `public_data/07_lanner_wot_full_visualization.ipynb` |
| Fig. 3B | Annotation and visualization of Day 4 LARRY-barcoded HCEBs | `LARRY/01_hceb_larry_annotation_visualization.ipynb` |
| Fig. 3C–D | Day 4 clone-sharing and permutation-based enrichment analysis | `LARRY/03_hceb_larry_basic_clonality.ipynb`; `LARRY/04_hceb_larry_clone_visualization.ipynb` |
| Fig. 3E–F | Day 4 endoderm subclustering and clone sharing with Trans ExM | `HCEB/05_endo_subtype.ipynb`; `LARRY/04_hceb_larry_clone_visualization.ipynb` |
| Fig. 3H–I | Day 7 LARRY annotation and clone-sharing analysis | `LARRY/01_hceb_larry_annotation_visualization.ipynb`; `LARRY/03_hceb_larry_basic_clonality.ipynb`; `LARRY/04_hceb_larry_clone_visualization.ipynb` |
| Fig. 3J–K | Day 12 LARRY annotation and clone-sharing analysis | `LARRY/01_hceb_larry_annotation_visualization.ipynb`; `LARRY/03_hceb_larry_basic_clonality.ipynb`; `LARRY/04_hceb_larry_clone_visualization.ipynb` |
| Extended Data Fig. 3A–C | Day 4 LARRY clone-size QC and epiblast clone-sharing analysis | `LARRY/03_hceb_larry_basic_clonality.ipynb`; `LARRY/04_hceb_larry_clone_visualization.ipynb` |
| Extended Data Fig. 3E–F | Day 7 clone-size QC and clone-sharing enrichment | `LARRY/03_hceb_larry_basic_clonality.ipynb`; `LARRY/04_hceb_larry_clone_visualization.ipynb` |
| Extended Data Fig. 3G | Day 12 clone-size QC | `LARRY/03_hceb_larry_basic_clonality.ipynb` |
| Extended Data Fig. 6B–D | WT versus CDX2-KO peri-gastruloid single-cell integration, annotation, and functional/module-score analyses | `PeriGas/02_integration.ipynb`; `PeriGas/03_annotation.ipynb`; `PeriGas/04_wt_ko.ipynb`; `PeriGas/05_function.ipynb` |
| Extended Data Fig. 6E | Reanalysis of human blastocyst bulk RNA-seq data | `public_data/03_icm_bulk.R` |
| Extended Data Fig. 8B–C | WT versus CDX2-KO differential expression and visualization of LMO2 | `PeriGas/04_wt_ko.ipynb`; `PeriGas/05_function.ipynb` |
| Fig. 6C | CRISPRa sgRNA preprocessing, MAGeCK enrichment analysis, and candidate ranking | `CRISPRa/01_CRISPRa_preprocess.sh`; `CRISPRa/02_CRISPRa_downstream_analysis.R` |


## Contact

Yiming Chao, Hongji Li, Rio Sugimura

Email: chym@connect.hku.hk, troyli990601@connect.hku.hk, rios@hku.hk
