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
│   ├── 01_lcxyw_construction.ipynb
│   ├── 02_lcxyw_hceb_exm_mapping.ipynb
│   └── 03_cs6_hceb_exm_mapping.ipynb
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
│   ├── 06_tvae_stage_wt.ipynb
│   └── 07_lanner_wot_full_visualization.ipynb
│
└── CRISPRa/
    ├── 01_CRISPRa_preprocess.sh
    ├── 02_CRISPRa_downstream_analysis.R
    └── 03_CRISPRa_screen_library_summary.xls
```

Contact: Yiming Chao, Hongji Li, Rio Sugimura

Email: chym@connect.hku.hk, troyli990601@connect.hku.hk, rios@hku.hk
