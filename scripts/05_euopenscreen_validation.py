"""Step 05 — EU OpenScreen validation of the ChEMBL pathogen models.

Evaluates how well our step-04 Ersilia model predictions
(``output/04_ersilia_predictions/euopenscreen/{eosid}.csv``, headline column ``consensus_score``)
rank the actives in the EU OpenScreen screening library — compounds the models never saw in
training. Analyses:

  1. Own primary assay: each of the 7 organisms with an EU OpenScreen primary assay, scored by
     its own model (AUROC + ROC curves + enrichment factor).
  1b. Own secondary assay: the same models on the merged secondary (confirmatory / dose-response)
     assays, for a primary-vs-secondary AUROC comparison.
  3. Shared vs exclusive hits: AUROC on exclusive actives (hit in 1 of the 7 primary assays, i.e.
     organism-specific) vs shared/non-exclusive actives (hit in >=2, i.e. pan-active), each against
     the same primary inactives.
  4. Cross-organism: model x EU OpenScreen assay AUROC matrix (off-diagonal = a model predicting
     a DIFFERENT organism's data) + per-model specificity index.
  4b. Active-set overlap: pairwise Jaccard between the 7 assays' active sets (label-only) — the
     backdrop for reading the cross-organism AUROCs.
  4c. Hit promiscuity (label-only): how many actives are hits in 1, 2, ... 7 pathogens, plus the
     per-compound table naming the promiscuous (pan-active) ones.
  4d. Summed consensus score across the 7 models, as a boxplot over four hit classes — inactive in
     every assay it was tested in / hit in exactly 1 pathogen / narrow-spectrum (2-3 pathogens) /
     broad-spectrum (>3). Primary assays only, training-set compounds included (raw, no leakage
     filtering). Plus the same data as MAXIMUM score over the 7 models, inactive vs active
     (a hit in >=1 pathogen), on both the raw score scale and after converting each model's score
     to its within-model library percentile (the raw scores are not calibrated across models).
  4e. Own-model rank for exclusive hits: for each hit in exactly 1 pathogen, the rank its own
     pathogen's model gives it among all 7 models (1 = its own pathogen scores it highest).
     Two panels — raw consensus scores, and scores converted to within-model percentiles first
     (the raw scores are not calibrated across models).
  + Per-submodel breakdown (every output column, not just consensus) and a leakage report.
  + eos3dys sub-analysis (subfolder ``eos3dys_validation/``): the CoAdd model eos3dys, whose many
    endpoints are scored against every EU OpenScreen assay (endpoint x assay heatmap + a
    same-vs-different-organism swarm). This tests the CoAdd model's generalization to EU OpenScreen
    — the opposite direction from the ChEMBL-models-on-CoAdd step (src/eval_coadd.py). Plus two
    own-assay analyses mirroring the ChEMBL ones, over the 6 organisms with BOTH an EU OpenScreen
    assay and an eos3dys endpoint (efaecium has no eos3dys endpoint): CoAdd-training overlap with
    the EU OpenScreen library and its actives, and exclusive-vs-shared hit AUROC per endpoint.
    Labels are always EU OpenScreen; CoAdd is only the dedup source.

Every metric (AUROC, AUPRC, BEDROC, EF@1%, EF@5%) is reported both ``raw`` and InChIKey-
``dedup`` (leakage-filtered against the ChEMBL training sets; for eos3dys, against its own CoAdd
training compounds). Missing step-04 model files are skipped with a logged message. Reads local
step-04 outputs and staged labels only; writes nothing outside output/05_euopenscreen_validation.

The ChEMBL-models-on-CoAdd validation is a separate concern and lives in its own step (see
src/eval_coadd.py + src/plots_coadd.py).

Output
------
Everything lands under ``output/05_euopenscreen_validation/``, filed by LEAKAGE STATUS so no number
is ambiguous about whether the models had seen the compound in training. Each of the three dirs
carries its own ``png/``, ``pdf/`` and ``figure_cells.json``.

    <top level>              — no leakage dimension: label-only tables + the leakage audit
        05_leakage_report.csv, 05_active_overlap.csv, 05_hit_promiscuity.csv,
        05_promiscuous_hits.csv
        figures: euos_overlap, active_overlap_jaccard, hit_promiscuity

    full/                    — training-set compounds KEPT
        05_euopenscreen_auroc.csv, 05_euopenscreen_secondary_auroc.csv, 05_euopenscreen_roc.csv,
        05_hit_exclusivity.csv, 05_cross_organism_euos.csv          (set=raw rows)
        05_consensus_sum_{boxstats,actives}.csv,
        05_consensus_max_{boxstats,actives}.csv,
        05_consensus_max_percentile_{boxstats,actives}.csv,
        05_exclusive_hit_model_rank{,_compounds}.csv
        figures: consensus_sum_by_hit_class, consensus_max_by_activity,
                 consensus_max_percentile_by_activity, exclusive_hit_model_rank_{raw,percentile}

    deduplicated/            — training-set compounds REMOVED
        05_euopenscreen_auroc.csv, 05_euopenscreen_secondary_auroc.csv, 05_euopenscreen_roc.csv,
        05_hit_exclusivity.csv, 05_cross_organism_euos.csv          (set=dedup rows)
        05_specificity_index.csv,
        05_consensus_max_percentile_dedup_{boxstats,actives}.csv,
        05_exclusive_hit_model_rank_dedup{,_compounds}.csv
        figures: euos_roc_grid, euos_shared_enrichment, primary_vs_secondary_auroc,
                 hit_exclusivity_auroc, cross_organism_heatmap, specificity_index,
                 submodel_auroc_summary, consensus_max_percentile_by_activity_dedup,
                 exclusive_hit_model_rank_percentile_dedup

The five long-form metric tables appear in BOTH subfolders, filtered to that folder's ``set``.
The AUROC-family figures sit under ``deduplicated/`` because they plot the leakage-filtered values;
their ``full/`` counterparts exist as data, not as figures.

``individual_performance/`` and ``eos3dys_validation/`` are separate analysis families and keep
their own layout (their CSVs carry a ``set`` column with both variants):

    output/05_euopenscreen_validation/individual_performance/05_submodel_auroc.csv
    output/05_euopenscreen_validation/individual_performance/{code}_submodel_corr.csv
    output/05_euopenscreen_validation/individual_performance/{png,pdf}/<panel>.{png,pdf}
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_euos_auroc.csv
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_overlap_report.csv
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_hit_exclusivity.csv
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_exclusive_rank.csv
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_exclusive_rank_compounds.csv
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_roc.csv
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_consensus_max_boxstats.csv
    output/05_euopenscreen_validation/eos3dys_validation/eos3dys_consensus_max_actives.csv
    output/05_euopenscreen_validation/eos3dys_validation/{png,pdf}/<panel>.{png,pdf}
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from eval_euopenscreen import run_all  # noqa: E402
from plots_euopenscreen import (  # noqa: E402
    save_euopenscreen_figures,
    save_individual_performance_figures,
)
from eval_eos3dys import run_all as run_eos3dys_all  # noqa: E402
from plots_eos3dys import save_eos3dys_figures  # noqa: E402

# Predictions from step 04 (eosvc-tracked, read locally).
pred_dir = os.path.join(root, "..", "output", "04_ersilia_predictions")
# Ground-truth labels staged by 00_download_data.py.
euos_root = os.path.join(root, "..", "data", "raw", "euopenscreen_data")
coadd_root = os.path.join(root, "..", "data", "raw", "coadd_data")
# Sibling ChEMBL models repo (leakage/dedup training sets; degrades gracefully if absent).
models_root = os.path.join(root, "..", "..", "chembl-antimicrobial-models")
config_path = os.path.join(root, "..", "config", "pathogens_of_interest.csv")
output_dir = os.path.join(root, "..", "output", "05_euopenscreen_validation")
os.makedirs(output_dir, exist_ok=True)

run_all(
    pred_dir=pred_dir,
    euos_root=euos_root,
    models_root=models_root,
    config_path=config_path,
    output_dir=output_dir,
)
save_euopenscreen_figures(output_dir)
save_individual_performance_figures(os.path.join(output_dir, "individual_performance"))

# CoAdd model (eos3dys) validated on EU OpenScreen — its own subfolder (dedup removes each
# endpoint's CoAdd-training compounds). Distinct from the ChEMBL-models-on-CoAdd step.
eos3dys_dir = os.path.join(output_dir, "eos3dys_validation")
os.makedirs(eos3dys_dir, exist_ok=True)
run_eos3dys_all(pred_dir, euos_root, coadd_root, config_path, eos3dys_dir)
save_eos3dys_figures(eos3dys_dir)

print(f"\nDone → {output_dir}")
