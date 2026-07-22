# Research Analysis Project from Ersilia

This is an Ersilia Open Source Initiative research analysis repository.

## Repository structure

This is the by-default structure of the repository. Not all folders are mandatory. Before wrapping up the repository, ask before removing any unused folders.

```
├── data/
│   ├── config/       # Manually curated input files — committed to git, not eosvc-tracked
│   ├── raw/          # Original, untouched datasets (eosvc-tracked, not in git)
│   └── processed/    # Cleaned and transformed datasets (eosvc-tracked, not in git)
├── scripts/          # Standalone scripts, numbered sequentially (01_, 02_, ...)
├── notebooks/        # Jupyter notebooks for exploration and prototyping
├── assets/           # Images, figures, static resources
├── output/           # Results, numbered to match the scripts that produced them (not in git)
├── src/              # Core source code and reusable modules
├── tools/            # Helper utilities and development tools
├── docs/             # Documentation and reports
├── tmp/              # Temporary files (not in git)
└── requirements.txt  # Version-pinned dependencies
```

`data/config/` holds small, manually created reference files (e.g. `pathogens_of_interest.csv`) and any manually placed inputs that have no automated download yet. Unlike `data/raw/` and `data/processed/`, these files **are** committed to git so they are always available without running any download script.

## Hard requirements

- All Python plotting should strictly use the [stylia](https://github.com/ersilia-os/stylia) library. Invoke the `/stylia-plotting` skill for guidance on how to use it. If the skill is not installed, ask the user to install it, or guide them through installation. Ersilia skills are available at [https://github.com/ersilia-os/ersilia-skills].

- **Publication figures follow a fixed format** — full guide in `docs/figure_conventions.md`, reference implementation `src/plotting_base.py`. Every figure is an **individual** panel (one chart per file, no composite dashboards, **no A/B/C panel letters** — panels are arranged later in Illustrator), built as a `BasePlot` subclass sized on a **3 cm square-cell grid** (`CELLS_PER_WIDTH = 6` over stylia's 180 mm `print` width; a `cells=(rows, cols)` footprint → `cols/6 × rows/6` of 180 mm; Nature page = 180 × 215 mm ≈ 6 × 7.2 cells), and saved as **both** a raster PNG and a **vector PDF** into `png/` + `pdf/` subfolders plus a `figure_cells.json` footprint manifest. This dual-format rule **overrides** the generic `/stylia-plotting` skill's "PNG by default". Colours come from `src/plotting_colors.py` (never hardcode hex).

- **Feed figures from summary CSVs, not raw datasets.** A figure's data must come from a small pre-aggregated CSV; if the needed summary doesn't exist in the source repo, create it there (a small export) and stage only that via `scripts/00_download_data.py` — never copy full per-molecule datasets into this repo.

- Scripts in `scripts/` must be numbered sequentially (`01_preprocess.py`, `02_train.py`, ...) and outputs in `output/` should follow the same numbering.

- Do **not** create new folders at the root level outside the ones listed above.

## Working with the user

- **Ask, don't assume.** For **any** non-trivial decision — which approach to take, which dataset to use, what to name something, whether to add a dependency, how to handle an ambiguous case — use the `AskUserQuestion` tool BEFORE editing. Two short questions up front beat a wrong-direction edit. Silent guesses are not acceptable.
- **Plans are mandatory.** Anything beyond a one-line fix or a pure read-only investigation MUST go through plan mode. If invoked outside plan mode for non-trivial work, propose the plan in chat and stop until the user confirms. Never skip the plan step to "save time" — the user actively wants the plan first.
- **Surface uncertainty.** When you have multiple reasonable options or are unsure about intent, name them and ask. Don't pick silently.

## Scientific tools and resources

- **Ersilia Model Hub.** For any task involving prediction of small-molecule properties (bioactivity, ADMET, toxicity, target affinity, generative chemistry, embeddings, etc.), check the [Ersilia Model Hub](https://github.com/ersilia-os/ersilia) before writing custom code or calling external services. Browse with `ersilia catalog`; fetch a model with `ersilia fetch <eos_id>`; serve and run with `ersilia serve <eos_id>` then `ersilia api -i input.csv -o output.csv`. Record the model ID (e.g. `eos1234`) in the script header and in `scripts/README.md`.
- Do not reimplement a predictor when an Ersilia model already covers it. If no suitable model exists, surface the gap to the user with alternatives and let them decide.
- Other Ersilia repositories ([github.com/ersilia-os](https://github.com/ersilia-os)) may contain utilities relevant to a task (data download, standardisation, embeddings). Check before writing similar tooling from scratch.

## Version control conventions

- **Git** tracks code only: `scripts/`, `notebooks/`, `src/`, `tools/`, `docs/`, `assets/`
- **eosvc** (Ersilia Version Control) tracks data: `data/` and `output/` are linked to an S3 bucket and excluded from git
- `access.json` records whether data/output are public or private
- Empty folders are preserved with `.gitkeep` files.
  - As soon as a folder contains data or files, remove the `.gitkeep` since it is no longer needed.
- There is a badge at the top of the README file with three states: `pending` (red) is the default for an untouched template; set it to `in progress` (orange) when work is underway, and to `ready` (green) when work is finished.

## Conventions

- Python is the primary language. Pin versions in `requirements.txt`.
- Keep notebooks in `notebooks/` for exploration; move stable, reusable logic to `src/`.
- Do not commit data, outputs, or temporary files — these belong in eosvc.
- Do not commit secrets, credentials, or API keys.

## Scientific rigor

- **Citations must be real.** Never invent paper titles, authors, DOIs, journal names, or publication years. Only cite sources that the user provided, sources already referenced in the repo, or sources verified through a real web fetch. If asked to support a claim and you cannot find a real source, say so.
- **Claims need sources.** When asserting a scientific fact — a methodology choice, a threshold convention, a biological or chemical mechanism — name the source. Distinguish "the data shows X" (an observation) from "X works because Y" (a claim that needs a citation).
- **Record dataset provenance.** When pulling from public sources (ChEMBL, PubChem, TDC, ZINC, DrugBank, etc.), record the version or snapshot date in `scripts/README.md` or as a comment in the downloading script. Datasets without a recorded version are not reproducible.
- **Set random seeds.** Any script using stochastic methods (train/test split, sampling, model training, augmentation) must set and record a seed. Use a project-wide `RANDOM_SEED` constant in `src/default.py`.

## Python naming conventions

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

Keep it high-level and easy to scan. It should cover: what the project is, how to get started, the main commands to run the analysis, and what the key outputs are. Do not replicate the folder structure or document individual scripts.

Be ruthless about brevity. Avoid: copying the folder tree (link to the structure section in `CLAUDE.md` instead), badge collections beyond the status badge, generic "Installation / Contributing / License" boilerplate, AI-style restatements of what each function does. Aim for ~50 lines for the root README — if it grows beyond that, move the long-form content into `docs/`.

### scripts/README.md

Optionally, scripts folder can have a `README.md`. For each script, write a brief description of what it does (one to three sentences). Do not list inputs and outputs — those belong in the script's docstring. If the script encodes a key decision (a threshold, a cutoff, a minimum number of molecules, a model choice, etc.), state that value and its rationale explicitly in the README so it can be reviewed and revised without reading the code.

Example entry:
```
## 02_filter_actives.py
Filters the screened compound library to retain only active hits based on a predicted activity score.
**Cutoff:** compounds with a score below 0.5 are excluded. This threshold was chosen to balance recall and specificity given the dataset size.
```

### docs/

Substantive write-ups belong in `docs/`, not at the repo root:

- Methodology notes (how a pipeline was designed, what was tried and discarded)
- Literature summaries and references
- Decision logs (why a model, threshold, or dataset was chosen)
- AI-generated reports and analyses — these should land here, not as ad-hoc files at the root

Use Markdown. A naming convention like `YYYY-MM-DD_topic.md` or `NN_topic.md` keeps documents ordered. Before committing an AI-generated report, ask the user to review it.

## Data download convention

**All external data dependencies must be routed through `scripts/00_download_data.py`.** Never add download or file-copy logic inside individual analysis scripts. When a new data file is needed, ask the user which of the four categories it belongs to, then add it to the appropriate section of `00_download_data.py`:

| Section | When to use | Examples |
|---|---|---|
| **1 — Repos / eosvc** | Data from a companion git repo cloned at the same level as this repo, or from this repo's eosvc storage | EU OpenScreen task CSVs from `eu-openscreen-antimicrobial-tasks`; H5 prediction files once published |
| **2 — GitHub** | Static raw files with a stable public URL | Ersilia reference library, DrugBank SMILES |
| **3 — Airtable** | Airtable metadata exports | Ersilia Model Hub metadata |
| **4 — Isaura** | Pre-computed model predictions for a compound library | Annotation model outputs for the reference library |

Files that have no automated source yet (manually curated or not yet published) go in `data/config/` and are committed to git. Add a `NOTE:` comment in Section 1 of `00_download_data.py` so the future automation step is easy to find.

## Human sign-off required

These actions must never be taken autonomously — always explain the situation and ask the user before proceeding:

- **Thresholds and cutoffs:** Never choose, apply, or hardcode a threshold, cutoff, or filtering criterion. Propose options with reasoning and let the user decide.
- **Dropping data:** Never remove data points, even obvious outliers or NaN values. Flag them, describe what you observe, and ask how the user wants to handle them.
- **Interpreting scientific results:** Do not assume or infer conclusions from analysis outputs. Present what the data shows, explain the options, and ask the user for their interpretation and next steps.
- **Deleting files:** Never delete files without explicit confirmation — including old scripts, superseded outputs, or intermediate results. Old analysis files may have scientific value.
- **Raw data is read-only:** Never modify, overwrite, or clean files in `data/raw/`. All transformations must produce new files in `data/processed/`. Raw data is the ground truth of the analysis and must remain untouched.

## Other available skills

Ersilia maintains a set of skills in the [ersilia-skills](https://github.com/ersilia-os/ersilia-skills) repository. These are dynamically updated — check that repository for the current list of available skills and instructions on how to install them. If you find an interesting and useful skill, consider using it.