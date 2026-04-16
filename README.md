# Blood originates in hypoblasts during embryonic development
Blood is essential for oxygen supply throughout life. The emergence of blood in the human embryo remains poorly understood. Our study leverages multiple stem cell embryo models and advanced lineage barcoding to unveil that hypoblast, originally regarded as forming the yolk sac wall, is heterogeneous and contributes to CDX2+ extraembryonic mesoderm, followed by hemoglobin+ cells as the first blood cells. CDX2 marks the hypoblast-to-hemoglobin+ cell trajectory that functionally sustains oxygen levels in embryo models. These hemoglobin+ cells molecularly and functionally resemble phagocytes. We show that the erythro-core regulatory network is poised in hypoblasts, and its boost endows erythropoiesis to both hypoblasts and phagocytes. Hypoblasts are the origin of the first blood in humans and non-human primates, providing a conceptual framework that earlier blood generation than expected fills the gap in the establishment of circulation. Further, the hypoblast is a place where primates may repurpose the phagocyte program to carry oxygen throughout embryos.

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

## Analysis scripts in this repo:

## Analysis Scripts Structure

```
scripts/
│
├── public_data/
│   ├── 01_human_embryo_reanalysis.R
│   ├── 02_mouse_embryo_reanalysis.R
│   └── 03_icm_bulk.R
│
├── HCEB/
│   ├── 01_qc.ipynb
│   ├── 02_integration.ipynb
│   ├── 03_annotation.ipynb
│   ├── 04_blood_trajectory.ipynb
│   ├── 05_endo_subtype.ipynb
│   └── 06_tvae_stage_d4_10.ipynb
│
├── LARRY/
│   ├── 01_qc.ipynb
│   ├── 02_annotation.ipynb
│   ├── 03_barcode_preprocess.ipynb
│   ├── 04_clone_analysis.ipynb
│   └── 05_lineage.py
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

Contact: Yiming Chao, Hongji Li, Rio Sugimura

Email: chym@connect.hku.hk, troyli990601@connect.hku.hk, rios@hku.hk
