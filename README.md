# eos-analysis-template

This repository provides a structured template for setting up new research analysis in Ersilia.

## Background

Replace this paragraph with a short description of the project. This description should explain the background or context of the project, specifying collaborators.

## Tracking details

The project is tracked by Git (mainly for code) and DVC (mainly for data):

* Tracked by Git and linked to a Github repository: only src, scripts and notebooks.
* Tracked by DVC and linked to a Google Drive folder inside "Projects/<<Repository name>>".

## Repository structure

This repository is organized as follows:

```
eos-analysis-template/
│
├── LICENSE
├── README.md
├── .gitignore
├── install.sh
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

This repository may contain data and outputs that are not stored in GitHub. You can use [evc](https://github.com/ersilia-os/evc) to download these files or, otherwise, simply download them and place them in the current folder following [this link](https://example.com).

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
