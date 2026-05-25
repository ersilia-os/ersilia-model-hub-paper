import matplotlib.colors as mcolors
import numpy as np
import stylia as st
from matplotlib.patches import Rectangle
from sklearn.metrics import roc_curve, auc


def abbrev(name):
    """Abbreviate genus: 'Mycobacterium tuberculosis' → 'M. tuberculosis'."""
    parts = name.split()
    return f"{parts[0][0]}. {' '.join(parts[1:])}" if len(parts) > 1 else name


def plot_roc_single(ax, y_true, y_pred, title, color=None):
    """Single ROC curve with AUC and class-count annotation.

    Parameters
    ----------
    ax     : matplotlib Axes
    y_true : array-like of int
    y_pred : array-like of float
    title  : str
    color  : optional color override (default nc.plum)

    Returns AUC (float) or None if only one class present.
    """
    nc = st.NamedColors()
    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    c = color if color is not None else nc.turquoise

    if len(np.unique(y_true)) < 2:
        ax.text(0.5, 0.5, "No positive\nexamples", ha="center", va="center",
                transform=ax.transAxes)
        st.label(ax, xlabel="", ylabel="", title=title)
        return None

    fpr, tpr, _ = roc_curve(y_true, y_pred)
    roc_auc = auc(fpr, tpr)
    n_pos = int(y_true.sum())
    n_neg = int(len(y_true) - n_pos)

    ax.plot(fpr, tpr, color=c, label=f"AUC={roc_auc:.2f}\n({n_pos}+ {n_neg}−)")
    ax.fill_between(fpr, tpr, alpha=0.15, color=c)
    ax.plot([0, 1], [0, 1], "--", color=nc.silver)
    st.label(ax, xlabel="FPR", ylabel="TPR", title=title)
    ax.legend(fontsize=6, loc="lower right")
    return roc_auc


def plot_auroc_dotplot(ax, auroc_df, title="Prediction performance — all assays"):
    """Horizontal dot plot of per-pathogen AUROC across all assay features.

    One dot per (pathogen, feature) pair, coloured by feature. No legend.
    Pathogens sorted ascending by consensus_score AUROC (or mean if absent).
    A dashed reference line marks 0.5.

    Parameters
    ----------
    ax       : matplotlib Axes
    auroc_df : DataFrame with columns pathogen, feature, auroc
    title    : str, plot title
    """
    features = list(auroc_df["feature"].unique())
    pal = st.CategoricalPalette("ersilia")
    colors = pal.get(len(features))
    feature_color = dict(zip(features, colors))

    if "consensus_score" in features:
        sort_series = (
            auroc_df[auroc_df["feature"] == "consensus_score"]
            .set_index("pathogen")["auroc"]
        )
    else:
        sort_series = auroc_df.groupby("pathogen")["auroc"].mean()
    pathogens_order = sort_series.sort_values(ascending=True).index.tolist()
    pathogen_idx = {p: i for i, p in enumerate(pathogens_order)}

    for feat in features:
        sub = auroc_df[auroc_df["feature"] == feat]
        ys = [pathogen_idx[p] for p in sub["pathogen"]]
        ax.scatter(sub["auroc"].values, ys, color=feature_color[feat])

    ax.axvline(0.5, color="gray", linestyle="--")
    ax.set_yticks(range(len(pathogens_order)))
    ax.set_yticklabels([abbrev(p) for p in pathogens_order])
    ax.set_xlim(0, 1)
    st.label(ax, xlabel="AUROC", ylabel="", title=title)



def plot_rank_boxplots(ax, models, title):
    """Horizontal paired boxplots of predicted ranks for actives vs inactives.

    One pair per model (actives = nc.turquoise, inactives = nc.crimson).
    Uses fold 0 predictions only. Models sorted ascending by mean_auroc (best at top).

    Parameters
    ----------
    ax     : matplotlib Axes
    models : list of dicts with keys: name, mean_auroc, rank_actives, rank_inactives
    title  : str — pathogen name used as panel title
    """
    nc = st.NamedColors()
    sorted_models = sorted(models, key=lambda x: x["mean_auroc"])
    n = len(sorted_models)

    for i, m in enumerate(sorted_models):
        y_pos = i * 2
        bp_a = ax.boxplot(
            m["rank_actives"], positions=[y_pos + 0.3], vert=False, widths=0.5,
            patch_artist=True, manage_ticks=False,
        )
        bp_i = ax.boxplot(
            m["rank_inactives"], positions=[y_pos - 0.3], vert=False, widths=0.5,
            patch_artist=True, manage_ticks=False,
        )
        for patch in bp_a["boxes"]:
            patch.set_facecolor(nc.turquoise)
        for patch in bp_i["boxes"]:
            patch.set_facecolor(nc.crimson)
        for element in ("whiskers", "caps", "medians", "fliers"):
            for line in bp_a[element]:
                line.set_color(nc.turquoise)
            for line in bp_i[element]:
                line.set_color(nc.crimson)

    ax.set_yticks([i * 2 for i in range(n)])
    ax.set_yticklabels([m["name"] for m in sorted_models])
    st.label(ax, xlabel="Predicted rank", ylabel="", title=title)


def plot_auroc_heatmap(ax, matrix_df, title="AUROC heatmap", highlight_cells=None):
    """Heatmap of AUROC values: rows = index labels, columns = task codes.

    Cells are coloured with a diverging colormap centred at 0.5 (random baseline).
    NaN cells (insufficient data) are shown in silver.
    Each valid cell is annotated with its AUROC value.

    Parameters
    ----------
    ax              : matplotlib Axes
    matrix_df       : DataFrame with row labels as index, task codes as columns, AUROC values
    title           : str, plot title
    highlight_cells : list of (row_idx, col_idx) tuples to outline with a red frame, or None
    """
    nc = st.NamedColors()
    data = matrix_df.values.astype(float)
    row_labels = [abbrev(p) for p in matrix_df.index]
    col_labels = list(matrix_df.columns)
    nrows, ncols = data.shape

    cm = st.DivergingColormap("crimson_cobalt")
    cm.fit(np.array([0.0, 1.0]))

    rgba = np.zeros((nrows, ncols, 4))
    for i in range(nrows):
        for j in range(ncols):
            v = data[i, j]
            if np.isnan(v):
                rgba[i, j] = mcolors.to_rgba(nc.silver)
            else:
                rgba[i, j] = cm.transform(np.array([v]))[0]

    ax.imshow(rgba, aspect="auto")

    if highlight_cells:
        for (ri, ci) in highlight_cells:
            rect = Rectangle(
                (ci - 0.5, ri - 0.5), 1, 1,
                linewidth=2, edgecolor="red", facecolor="none",
            )
            ax.add_patch(rect)

    for i in range(nrows):
        for j in range(ncols):
            v = data[i, j]
            if not np.isnan(v):
                text_color = "white" if (v > 0.75 or v < 0.25) else "black"
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=st.FONTSIZE_SMALL, color=text_color)

    ax.set_xticks(range(ncols))
    ax.set_xticklabels(col_labels, rotation=45, ha="right")
    ax.set_yticks(range(nrows))
    ax.set_yticklabels(row_labels)
    st.label(ax, xlabel="", ylabel="", title=title)


def plot_same_vs_diff_swarm(ax, same_aurocs, diff_aurocs, title="Same vs different pathogen AUROC"):
    """Jittered swarm + boxplot comparing same-pathogen vs different-pathogen AUROC.

    Parameters
    ----------
    ax          : matplotlib Axes
    same_aurocs : array-like of float, AUROC values for same-pathogen (endpoint, task) pairs
    diff_aurocs : array-like of float, AUROC values for different-pathogen pairs
    title       : str
    """
    nc = st.NamedColors()
    rng = np.random.default_rng(42)

    groups = [np.array(same_aurocs, dtype=float), np.array(diff_aurocs, dtype=float)]
    colors = [nc.turquoise, nc.crimson]
    x_labels = ["Same\npathogen", "Different\npathogen"]

    for i, (vals, color) in enumerate(zip(groups, colors)):
        vals = vals[~np.isnan(vals)]
        jitter = rng.uniform(-0.15, 0.15, size=len(vals))
        ax.scatter(np.full(len(vals), i) + jitter, vals, color=color)
        bp = ax.boxplot(
            vals, positions=[i], vert=True, widths=0.3,
            patch_artist=True, manage_ticks=False, showfliers=False,
        )
        for patch in bp["boxes"]:
            patch.set_facecolor("none")
            patch.set_edgecolor(color)
        for element in ("whiskers", "caps", "medians"):
            for line in bp[element]:
                line.set_color(color)

    ax.axhline(0.5, color=nc.silver, linestyle="--")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(x_labels)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0, 1)
    st.label(ax, xlabel="", ylabel="AUROC", title=title)


def plot_specificity_bars(ax, df, title):
    """Horizontal bar chart of pathogen model specificity index.

    Bars coloured turquoise (positive index) or crimson (negative).
    A dashed reference line marks 0.

    Parameters
    ----------
    ax    : matplotlib Axes
    df    : DataFrame with columns 'pathogen' and 'specificity_index'
            (rows with NaN specificity_index are ignored)
    title : str
    """
    nc = st.NamedColors()
    df = df.dropna(subset=["specificity_index"]).sort_values(
        "specificity_index", ascending=True
    )
    labels = [abbrev(p) for p in df["pathogen"]]
    values = df["specificity_index"].values
    colors = [nc.turquoise if v >= 0 else nc.crimson for v in values]

    ax.barh(range(len(df)), values, color=colors)
    ax.axvline(0, color=nc.silver, linestyle="--")
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels)
    st.label(
        ax,
        xlabel="Specificity index (same − mean cross-pathogen AUROC)",
        ylabel="",
        title=title,
    )
