# Ersilia Analysis Project

This is an Ersilia Open Source Initiative research analysis repository. It follows the [eos-analysis-template](https://github.com/ersilia-os/eos-analysis-template) structure.

## Organisation context

[Ersilia](https://ersilia.io) is a tech-nonprofit building open-source AI/ML tools for antimicrobial drug discovery, focused on the Global South. The main asset is the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia).

## Repository structure

```
├── data/
│   ├── raw/          # Original, untouched datasets (eosvc-tracked, not in git)
│   └── processed/    # Cleaned and transformed datasets (eosvc-tracked, not in git)
├── scripts/          # Standalone scripts, numbered sequentially (01_, 02_, ...)
├── notebooks/        # Jupyter notebooks for exploration and prototyping
├── assets/           # Images, figures, static resources
├── output/           # Results, numbered to match the scripts that produced them (not in git)
├── src/              # Core source code and reusable modules
├── tools/            # Helper utilities and development tools
├── docs/             # Documentation, reports, AI-generated files
├── tmp/              # Temporary files (not in git)
└── requirements.txt  # Version-pinned dependencies
```

## Version control conventions

- **Git** tracks code only: `scripts/`, `notebooks/`, `src/`, `tools/`, `docs/`, `assets/`
- **eosvc** (Ersilia Version Control) tracks data: `data/` and `output/` are linked to an S3 bucket and excluded from git
- `access.json` records whether data/output are public or private
- Empty folders are preserved with `.gitkeep` files

## Coding guidelines

- Python is the primary language. Pin versions in `requirements.txt`.
- Scripts in `scripts/` must be numbered sequentially (`01_preprocess.py`, `02_train.py`, ...) and outputs in `output/` should follow the same numbering.
- Keep notebooks in `notebooks/` for exploration; move stable, reusable logic to `src/`.
- Do not commit data, outputs, or temporary files — these belong in eosvc.
- Do not commit secrets, credentials, or API keys.

## Python naming conventions

- Variables inside a script use `snake_case` and are never capitalised.
- Project-wide constants (values reused across scripts) must be defined in `src/default.py` and named in `ALL_CAPS`.
- Scripts that import from `src/` must include this path setup at the top, before any `src` imports:

```python
import os
import sys
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
```

- Declare input and output folder paths as variables at the top of the script (module level, not inside functions) and ensure they exist with `os.makedirs(..., exist_ok=True)`. Do not create folders inside functions unless strictly necessary for that function's logic.

```python
data_dir = os.path.join(root, "..", "data", "processed")
output_dir = os.path.join(root, "..", "output")
os.makedirs(data_dir, exist_ok=True)
os.makedirs(output_dir, exist_ok=True)
```

## README guidelines

### Root README
Keep it high-level and easy to scan. It should cover: what the project is, how to get started, the main commands to run the analysis, and what the key outputs are. Do not replicate the folder structure or document individual scripts — that detail belongs elsewhere.

### scripts/README.md
Every scripts folder must have a `README.md`. For each script, write a brief description of what it does (one to three sentences). Do not list inputs and outputs — those belong in the script's docstring. If the script encodes a key decision (a threshold, a cutoff, a minimum number of molecules, a model choice, etc.), state that value and its rationale explicitly in the README so it can be reviewed and revised without reading the code.

Example entry:
```
## 02_filter_actives.py
Filters the screened compound library to retain only active hits based on a predicted activity score.
**Cutoff:** compounds with a score below 0.5 are excluded. This threshold was chosen to balance recall and specificity given the dataset size.
```

## Human sign-off required

These actions must never be taken autonomously — always explain the situation and ask the user before proceeding:

- **Thresholds and cutoffs:** Never choose, apply, or hardcode a threshold, cutoff, or filtering criterion. Propose options with reasoning and let the user decide.
- **Dropping data:** Never remove data points, even obvious outliers or NaN values. Flag them, describe what you observe, and ask how the user wants to handle them.
- **Interpreting scientific results:** Do not assume or infer conclusions from analysis outputs. Present what the data shows, explain the options, and ask the user for their interpretation and next steps.
- **Deleting files:** Never delete files without explicit confirmation — including old scripts, superseded outputs, or intermediate results. Old analysis files may have scientific value.
- **Raw data is read-only:** Never modify, overwrite, or clean files in `data/raw/`. All transformations must produce new files in `data/processed/`. Raw data is the ground truth of the analysis and must remain untouched.

## Plotting

Use the `stylia` package for all plots (Ersilia's matplotlib wrapper for publication-ready figures). Invoke the `/stylia-plotting` skill for guidance on how to use it.

## Available skills

Ersilia maintains a set of Claude Code skills in the [ersilia-skills](https://github.com/ersilia-os/ersilia-skills) repository. These are dynamically updated — check that repository for the current list of available skills and instructions on how to install them.
