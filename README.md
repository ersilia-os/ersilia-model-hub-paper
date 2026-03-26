# Ersilia Model Hub figures and analyses

This repository provides a structured template for setting up new research analysis in Ersilia.

## Background

Replace this paragraph with a short description of the project. This description should explain the background or context of the project, specifying collaborators.

## Tracking details

Data artifacts (`data/`, `output/`) are synced to S3 via [eosvc](https://github.com/ersilia-os/eosvc) and are not stored in Git. Access rules are defined in `access.json` at the repo root.

## Repository structure

This repository is organized as follows:

```
eos-analysis-template/
│
├── LICENSE
├── README.md
├── .gitignore
├── requirements.txt
│
├── data/
│   ├── raw/
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
- **.git/** → Git metadata (version control)  

---

📌 Empty folders are preserved with `.gitkeep` files so the structure remains consistent in Git.

---

## Project motivation and goal

Write a brief description about the scientific motivation and goal of the project. 

## 🚀 Getting Started

1. **Clone this repository**  
   ```bash
   git clone <your-repo-url>
   cd eos-analysis-template
  ```

## Using this repository

Data and outputs are not stored in Git. Use `eosvc` to sync them from S3.

**Install:**
```bash
pip install -r requirements.txt
```

**Set up credentials** (skip if accessing public data only):
```bash
eosvc config --access-key-id "..." --secret-access-key "..." --region "eu-central-2"
```

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
