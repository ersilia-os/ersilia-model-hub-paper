"""Step 06 — CoAdd validation of the public ChEMBL pathogen models.

Mirror of step 05 (EU OpenScreen) on the CoAdd dataset: the same step-04 predictions
(``output/04_ersilia_predictions/coadd/{eosid}.csv``, headline column ``consensus_score``) are
scored against CoAdd growth-inhibition / MIC labels. This first version stays narrow — only the
organisms that have both a ChEMBL model and a CoAdd reference strain, scored on that reference
strain at the two headline cutoffs:

  - Own-strain AUROC: each of the 8 organisms in ``COADD_REF_STRAINS`` (abaumannii, calbicans,
    ecoli, efaecium, kpneumoniae, paeruginosa, saureus, spneumoniae) on its reference strain, for
    both endpoints — single-point ``inhib_50`` and MIC ``mic_10``. efaecium and spneumoniae are
    MIC-only (inhib_50 skipped, logged).
  - Leakage report: overlap between each model's ChEMBL training set and its CoAdd eval strain.

Every metric is reported both ``raw`` and InChIKey-``dedup`` (leakage-filtered against the ChEMBL
training sets). The CoAdd model eos3dys is NOT evaluated here — it is validated the other way
round (on EU OpenScreen) in step 05's ``eos3dys_validation/``. Missing step-04 model files are
skipped with a logged message, so the step runs on partial data and re-runs cleanly once step 04
finishes. The richer multi-strain / multi-cutoff CoAdd matrix is deferred to a later version.

Output
------
    output/06_coadd_validation/06_coadd_auroc.csv
    output/06_coadd_validation/06_coadd_leakage_report.csv
    output/06_coadd_validation/png/<panel>.png
    output/06_coadd_validation/pdf/<panel>.pdf
    output/06_coadd_validation/figure_cells.json
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from eval_coadd import run_all  # noqa: E402
from plots_coadd import save_coadd_figures  # noqa: E402

# Predictions from step 04 (eosvc-tracked, read locally).
pred_dir = os.path.join(root, "..", "output", "04_ersilia_predictions")
# CoAdd ground-truth labels staged by 00_download_data.py.
coadd_root = os.path.join(root, "..", "data", "raw", "coadd_data")
# Sibling ChEMBL models repo (leakage/dedup training sets; degrades gracefully if absent).
models_root = os.path.join(root, "..", "..", "chembl-antimicrobial-models")
config_path = os.path.join(root, "..", "config", "pathogens_of_interest.csv")
output_dir = os.path.join(root, "..", "output", "06_coadd_validation")
os.makedirs(output_dir, exist_ok=True)

run_all(
    pred_dir=pred_dir,
    coadd_root=coadd_root,
    models_root=models_root,
    config_path=config_path,
    output_dir=output_dir,
)
save_coadd_figures(output_dir)
print(f"\nDone → {output_dir}")
