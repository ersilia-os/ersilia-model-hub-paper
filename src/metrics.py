"""Shared classification metrics for the external-validation step (05).

Ported from ``new-modelling/src/metrics.py``. Kept dependency-light (numpy +
sklearn only) so the eval/plot steps don't load any heavy modelling stack.

The active-ranking metrics matter here because the external screening libraries
(EU OpenScreen, CoAdd) have very low hit rates (~0.01-0.4% active), where AUROC
alone is a weak signal:
  - AUPRC        : precision-recall area, sensitive to class imbalance.
  - BEDROC(a=20) : early-recognition emphasis (~top 8%); Truchon & Bayly, 2007.
  - EF@1%/@5%    : enrichment factor = hit rate in the top k% over the base rate.
"""

import numpy as np
from sklearn.metrics import average_precision_score, roc_auc_score

BEDROC_ALPHA = 20.0  # early-recognition emphasis (~top 8%); documented choice
EF_FRACTIONS = (0.01, 0.05)  # enrichment-factor top-fractions reported


def bedroc(y_true, y_score, alpha=BEDROC_ALPHA):
    """Boltzmann-Enhanced Discrimination of ROC (Truchon & Bayly, 2007)."""
    y_true = np.asarray(y_true, dtype=float)
    y_score = np.asarray(y_score, dtype=float)
    n = len(y_true)
    na = int(y_true.sum())
    if na == 0 or na == n:
        return np.nan
    order = np.argsort(-y_score, kind="mergesort")
    y = y_true[order]
    ranks = np.where(y == 1)[0] + 1
    ra = na / n
    rie_num = np.sum(np.exp(-alpha * ranks / n))
    rie_den = ra * (1 - np.exp(-alpha)) / (np.exp(alpha / n) - 1)
    rie = rie_num / rie_den
    factor = ra * np.sinh(alpha / 2) / (np.cosh(alpha / 2) - np.cosh(alpha / 2 - alpha * ra))
    return rie * factor + 1.0 / (1.0 - np.exp(alpha * (1 - ra)))


def enrichment_factor(y_true, y_score, frac):
    """Enrichment factor in the top ``frac`` of ranked predictions."""
    y_true = np.asarray(y_true)
    n = len(y_true)
    na = int(y_true.sum())
    if na == 0:
        return np.nan
    k = max(1, int(round(n * frac)))
    top = np.argsort(-np.asarray(y_score), kind="mergesort")[:k]
    hits = int(np.asarray(y_true)[top].sum())
    return (hits / k) / (na / n)


def compute_metrics(y, score):
    """AUROC / AUPRC / BEDROC / EF + counts. NaN metrics when degenerate."""
    y = np.asarray(y, dtype=int)
    score = np.asarray(score, dtype=float)
    na = int(y.sum())
    out = {
        "n_eval":     len(y),
        "n_active":   na,
        "prevalence": round(na / len(y), 5) if len(y) else np.nan,
    }
    if na == 0 or na == len(y):
        out.update({"auroc": np.nan, "auprc": np.nan, "bedroc": np.nan,
                    "ef_1pct": np.nan, "ef_5pct": np.nan})
        return out
    out["auroc"]   = round(roc_auc_score(y, score), 4)
    out["auprc"]   = round(average_precision_score(y, score), 4)
    out["bedroc"]  = round(bedroc(y, score), 4)
    out["ef_1pct"] = round(enrichment_factor(y, score, 0.01), 3)
    out["ef_5pct"] = round(enrichment_factor(y, score, 0.05), 3)
    return out
