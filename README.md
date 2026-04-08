# Ersilia Model Hub figures and analyses

Code, data, and figures for the Ersilia Model Hub paper.

## Project motivation and goal

This repository contains the data, analyses, and figures used for the Ersilia Model Hub paper (INSERT LINK OF PUBLICATION). Several analyses are based on the [Ersilia reference set of compounds](https://github.com/ersilia-os/ersilia-model-hub-maintained-inputs), and precalculated model outputs are stored and retrieved via [Isaura](https://github.com/ersilia-os/isaura), Ersilia's precalculation store built on top of S3-compatible object storage.

## Tracking details

- **Git** (this GitHub repository): `src/`, `scripts/`, `notebooks/`, `assets/`, `docs/`, `tools/`
- **eosvc** (S3 storage, not tracked by Git): `data/`, `output/`

Access rules for [eosvc](https://github.com/ersilia-os/eosvc) are defined in `access.json` at the repo root.

## Repository structure

This repository is organized as follows:

```
ersilia-model-hub-paper/
│
├── LICENSE
├── README.md
├── .gitignore
├── requirements.txt
│
├── data/
│   ├── raw/
│   │   ├── compounds/
│   │   ├── isaura/
│   │   └── ersilia_metadata.csv
│   └── processed/
│
├── scripts/
├── notebooks/
├── assets/
├── output/
│   ├── results/
│   └── plots/
│
├── src/
├── tools/
├── docs/
├── tmp/
│
└── .git/
```

- **data/**
  - **raw/** → Original, untouched datasets
    - **compounds/** → Compound sets used as model inputs (e.g. Ersilia reference library)
    - **isaura/** → Precalculated model outputs downloaded from Isaura
  - **processed/** → Cleaned and transformed datasets
- **scripts/** → Standalone scripts for preprocessing or automation  
- **notebooks/** → Jupyter notebooks for exploration and prototyping  
- **assets/** → Images, figures, and other static resources  
- **output/**
  - **results/** → Numerical results, logs, or text outputs  
  - **plots/** → Visualizations and charts  
- **src/** → Core source code and reusable modules  
- **tools/** → Helper utilities and development tools  
- **docs/** → Project documentation and reports  
- **tmp/** → Temporary files or intermediate outputs

---

📌 Empty folders are preserved with `.gitkeep` files so the structure remains consistent in Git.

---

## 🚀 Getting Started

```bash
git clone https://github.com/ersilia-os/ersilia-model-hub-paper
cd ersilia-model-hub-paper
conda create -n ersiliapaper python=3.10 -y
conda activate ersiliapaper
pip install -r requirements.txt
```

## Utilities

- **`scripts/download_reference.py`** — downloads the Ersilia reference library and fetches precalculated Isaura outputs for a given model from the cloud. Output is saved to `data/raw/isaura/` following the naming convention `emh_paper_<model_id>_<version>.csv`:

```bash
python scripts/download_reference.py --model eos42ez --version v1
```

- **`scripts/download_ersilia_metadata.py`** — downloads Ersilia model metadata from the public Airtable base and saves it to `data/raw/ersilia_metadata.csv`:

```bash
python scripts/download_ersilia_metadata.py
```

## Using this repository

Data and outputs are not stored in Git. Use `eosvc` to sync them from S3.

**Set up credentials** (skip if accessing public data only):
```bash
eosvc config --access-key-id "..." --secret-access-key "..." --region "eu-central-2"
```
For further details see the [eosvc repository](https://github.com/ersilia-os/eosvc).

**Download data:**
```bash
eosvc download --path data/
eosvc download --path output/
```

**Upload data:**
```bash
eosvc upload --path data/
eosvc upload --path output/
```

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
