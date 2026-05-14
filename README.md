# Ersilia Model Hub figures and analyses

Code, data, and figures for the Ersilia Model Hub paper.

## Project motivation and goal

This repository contains the data, analyses, and figures used for the Ersilia Model Hub paper (INSERT LINK OF PUBLICATION). Several analyses are based on the [Ersilia reference set of compounds](https://github.com/ersilia-os/ersilia-model-hub-maintained-inputs), and precalculated model outputs are stored and retrieved via [Isaura](https://github.com/ersilia-os/isaura), Ersilia's precalculation store built on top of S3-compatible object storage.

## Tracking details

- **Git** (this GitHub repository): `src/`, `scripts/`, `notebooks/`, `assets/`, `docs/`, `tools/`
- **eosvc** (S3 storage, not tracked by Git): `data/`, `output/`

Access rules for [eosvc](https://github.com/ersilia-os/eosvc) are defined in `access.json` at the repo root.

## 🚀 Getting Started

Make sure you have all requirements installed in a conda environment, the latest code pulled from Git and the data and outputs syced from S3 with `eosvc`. For further details see the [eosvc repository](https://github.com/ersilia-os/eosvc).

**Download data:**
```bash
eosvc download --path data/
eosvc download --path output/
```

## About the Ersilia Open Source Initiative

The [Ersilia Open Source Initiative](https://ersilia.io) is a tech-nonprofit organization fueling sustainable research in the Global South. Ersilia's main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia), an open-source repository of AI/ML models for antimicrobial drug discovery.

![Ersilia Logo](assets/Ersilia_Brand.png)
