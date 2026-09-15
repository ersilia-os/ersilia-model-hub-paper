"""Step 15 figure — one pathogen's own endpoint x endpoint top-1000 shared-actives matrix.

Step 09 (:mod:`eval_correlations`) computes and caches the full 307x307 top-1000 Jaccard matrix
across every selected bioactivity endpoint, but only ever consumes it aggregated into same-pathogen
/ different-pathogen box plots. Nothing draws the raw endpoint x endpoint matrix for a single
pathogen. This module does exactly that: the SAME diagonal-blanking and continuous sequential-scale
conventions as step 10's overlap matrix (:mod:`plots_auroc_matrix`), reused rather than imported
(those helpers are tied to that module's multi-track chrome, which a single-organism matrix does not
need).

**Raw shared-actives COUNT (user-directed, 2026-09-03), not Jaccard.** The cell value is how many of
the two endpoints' top-1000 highest-scoring compounds are the SAME compounds — the same quantity
step 10's own overlap matrix draws — not the Jaccard ratio the module used to show. The colour scale
is therefore a plain linear ``Normalize(0, top_n)`` rather than a ``[0, 1]`` scale.

**A second figure, a chord (circos-style) diagram, is an ALTERNATIVE view of the exact same E. coli
matrix (user-directed, 2026-09-03)** — no new data, same ``sub`` DataFrame :func:`run_pathogen`
already builds. The 24 bioactivity endpoints sit on a ring; a ribbon between two endpoints has the
SAME colour meaning as a matrix cell (coloured by shared-actives count, on the identical
``Normalize(0, top_n)`` scale) and is WIDE where it leaves each endpoint, pinched toward the centre —
real circos-style ribbons, not a uniform-width stroke (see the pyCirclize note below for why that
distinction needed a library). The predictor node (``abx`` for E. coli) is deliberately NOT a 25th
ribbon node — it is drawn as a radial bar track around the ring instead, one bar per endpoint, so
"how much does this endpoint's top-1000 look like a known antibiotic" reads as a single glance around
the circle rather than one more crowded ribbon. **Nothing is filtered**: a zero-count pair is simply
not drawn, so an endpoint that turns out to correlate with nothing in this set is exactly as visible
as one drawn correctly — as a ring position with no ribbons.

**Two ring-ordering modes (user-directed, 2026-09-03)**, see :func:`_ordered_labels`: the FIRST
version placed endpoints in the matrix's model-grouped order, which — with ribbons drawn between
essentially every pair — made the figure look like an arbitrary mandala rather than showing
structure. ``order="cluster"`` fixes that with hierarchical clustering (SciPy's own optimal-leaf-
ordering) on the shared-actives counts, so strongly-linked endpoints sit next to each other and
strong ribbons become short local arcs. ``order="abx"`` (by convention, the pathogen's own
``predictor_family``) instead ranks the ring directly by overlap with the predictor node, descending
clockwise from 12 o'clock — a different question (a ranking) answered independently of ribbon
structure.

**pyCirclize, an explicit exception to this repo's "stylia only" plotting rule (``CLAUDE.md``),
scoped to the chord diagram alone (user-directed, 2026-09-03).** A matplotlib line stroke has
uniform width along its length, so no amount of ribbon-ordering ever produces the wide-at-the-rim,
pinched-at-the-centre shape a circos ribbon needs — that shape comes from each ribbon claiming a
real ARC SPAN on both endpoints (proportional to that link's share of the endpoint's total
connectivity), which is a solved layout problem (github.com/moshi4/pyCirclize,
``Circos.chord_diagram``), not a good use of from-scratch geometry code. Two things this needed,
verified against toy data before use:

- ``Circos.chord_diagram`` treats its input matrix as DIRECTED — fed the full symmetric shared-count
  matrix it draws every pair TWICE (double-counted sector sizes, duplicate ribbons). The lower
  triangle is zeroed before the matrix is handed over (see :func:`PathogenChordPlot`), using a rule
  fixed by each label's position in ``matrix.columns`` — independent of whichever `order` is
  active — so the same pair is never dropped or duplicated regardless of ring order.
- ``link_kws_handler(from_label, to_label) -> dict`` overrides one link's colour at a time (each of
  our sectors already IS a single endpoint, never grouped), which is what makes continuous
  magnitude-based ribbon colouring possible on top of pyCirclize's own per-sector ``cmap``.

``Circos.plotfig`` accepts an existing ``PolarAxes`` — :func:`_new_polar_axes` builds one sized by
the same cells-grid maths ``BasePlot`` already uses (``stylia.create_figure`` has no polar support),
then hands it to ``BasePlot.__init__`` via its normal "ax was supplied" branch, so
``BasePlot.save()`` (PNG + PDF via ``stylia.save_figure``) needs no changes at all.
"""

import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import stylia as st
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pycirclize import Circos

from default import ACTIVITY_BINARIZE_TOP_N, OVERLAP_MATRIX_SPECTRUM
from plotting_base import CELLS_PER_WIDTH, BasePlot
from plotting_colors import INK, REFERENCE_LINE, distinct_colors
from plotting_utils import heatmap, merge_figure_cells, spectrum_cmap

#: A cell for a pathogen with only 24 endpoints at 180 mm width is still legible with printed
#: values; a full square page footprint keeps individual cells close to step 10's ~6 mm density.
CELLS = (6, 6)

TRACK_SIZE = "2.2%"
TRACK_PAD = 0.035

#: Extra pad (points) for the row (y) tick labels. The left-hand model track is inserted into
#: exactly the space matplotlib reserves for y tick labels, so without this the band is drawn over
#: the end of the label text — the same fix step 10's AurocMatrixPlot applies for its own left
#: track(s), sized here for a single band rather than step 10's two.
Y_LABEL_PAD = 20


def _dark_threshold(cmap, vmax, target_luminance=0.5, resolution=256):
    """The value at which ``cmap`` over ``[0, vmax]`` first drops below ``target_luminance``.

    Computed from the colormap itself (as :func:`plots_auroc_matrix._dark_threshold` does) rather
    than guessed as a fraction of ``vmax``, so annotation text stays legible if the palette changes.
    """
    t = np.linspace(0.0, 1.0, resolution)
    rgba = np.array([cmap(v) for v in t])
    luminance = 0.2126 * rgba[:, 0] + 0.7152 * rgba[:, 1] + 0.0722 * rgba[:, 2]
    below = np.flatnonzero(luminance < target_luminance)
    frac = t[below[0]] if len(below) else 1.0
    return frac * vmax


def _blank_diagonal(matrix):
    """Copy of ``matrix`` with its diagonal set to NaN — self-overlap is ``top_n`` by construction,
    a property of the axes rather than a measurement (same rationale as step 10's matrices)."""
    out = matrix.copy()
    for i in range(len(out)):
        out.iloc[i, i] = np.nan
    return out


def _mark_diagonal(ax, n, color=None):
    """Dashed outline on each blanked diagonal cell, so its position stays visible without colour."""
    color = INK if color is None else color
    for i in range(n):
        ax.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, edgecolor=color,
                               linestyle=(0, (1.5, 1.5)), linewidth=0.5, zorder=5))


def _model_color_map(model_id):
    """``{model: colour}``, in stable first-appearance order.

    Shared by every plot in this module (the matrix's left-hand track, the chord diagram's node
    dots) so the SAME model always gets the SAME colour across both E. coli figures — a reader who
    has seen one legend recognizes the other without re-reading it.
    """
    models = list(dict.fromkeys(model_id))  # stable order of first appearance
    return dict(zip(models, distinct_colors(len(models))))


class PathogenEndpointMatrixPlot(BasePlot):
    """One pathogen's own endpoint x endpoint top-1000 shared-actives matrix, with a model-colour
    track.

    ``matrix`` is the pathogen's own N x N submatrix of RAW shared-actives counts (index/columns
    already the readable ``"{model_id}:{column_name}"`` keys). ``model_id`` is a same-length
    sequence giving each row's source model, used only for the left-hand colour track and its
    legend. ``top_n`` sets the colour scale's ceiling and the colourbar label — it must match the
    cutoff the matrix was actually computed at.
    """

    def __init__(self, matrix, model_id, ax=None, cells=CELLS, top_n=ACTIVITY_BINARIZE_TOP_N,
                 name=None):
        super().__init__(ax=ax, cells=cells)
        self.name = name or "15_endpoint_overlap_matrix"
        if not matrix.size:
            self._unavailable()
            return

        shown = _blank_diagonal(matrix)
        cmap, norm = spectrum_cmap(OVERLAP_MATRIX_SPECTRUM), mcolors.Normalize(vmin=0, vmax=top_n)
        dark_at = _dark_threshold(cmap, top_n)

        heatmap(self.ax, shown, cmap=cmap, norm=norm, annotate=True, value_fmt="{:.0f}",
                text_light_when=lambda v: v >= dark_at, nan_color=REFERENCE_LINE,
                x_rotation=90, annot_fontsize=st.FONTSIZE_SMALL, aspect="equal")
        _mark_diagonal(self.ax, len(matrix))
        self.ax.tick_params(labelsize=st.FONTSIZE_SMALL, length=0)
        self.ax.tick_params(axis="y", pad=Y_LABEL_PAD)
        self.ax.grid(False)
        self.ax.set_axisbelow(False)
        self.ax.set_xlabel("")
        self.ax.set_ylabel("")

        model_colors = _model_color_map(model_id)
        models = list(model_colors)

        divider = make_axes_locatable(self.ax)
        self._model_track(divider, model_id, model_colors)
        self._colorbar(cmap, norm, divider, top_n)
        self.legend(model_colors, loc="upper center", bbox_to_anchor=(0.5, -0.32),
                   ncol=min(len(models), 4))

    def _model_track(self, divider, model_id, model_colors):
        """Vertical colour band on the left, one swatch per row's source model.

        Shares y with the main axes (``sharey``) so the band stays cell-aligned regardless of how
        the divider has already shrunk ``self.ax`` — the same guarantee step 10's tracks rely on.
        """
        ax = divider.append_axes("left", size=TRACK_SIZE, pad=TRACK_PAD, sharey=self.ax)
        rgba = np.array([mcolors.to_rgba(model_colors[m]) for m in model_id])
        n = len(model_id)
        ax.imshow(rgba.reshape(-1, 1, 4), aspect="auto", interpolation="none",
                  extent=(-0.5, 0.5, n - 0.5, -0.5))
        ax.xaxis.set_visible(False)
        ax.yaxis.set_visible(False)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(False)
        return ax

    def _colorbar(self, cmap, norm, divider, top_n):
        from matplotlib.cm import ScalarMappable

        cax = divider.append_axes("right", size="1.2%", pad=0.12)
        cb = self.fig.colorbar(ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                               orientation="vertical")
        cb.set_label(f"Shared actives (of {top_n})", fontsize=st.FONTSIZE_SMALL)
        cb.ax.tick_params(labelsize=st.FONTSIZE_SMALL)


def save_pathogen_endpoint_matrix_figure(output_dir, matrix, model_id, name,
                                         top_n=ACTIVITY_BINARIZE_TOP_N):
    """Build the matrix figure and merge its footprint into ``figure_cells.json``.

    Routes through :func:`plotting_utils.merge_figure_cells` (not a direct write) so a future second
    pathogen figure drawn into this same output dir accumulates rather than truncates the manifest —
    the same convention every other multi-figure step in this repo follows.
    """
    plot = PathogenEndpointMatrixPlot(matrix, model_id, name=name, top_n=top_n)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print(f"  [skip figure] {plot.name}: empty matrix")
    return merge_figure_cells(output_dir, footprints)


# --------------------------------------------------------------------------- #
# Chord (circos-style) diagram — an alternative view of the same matrix       #
# --------------------------------------------------------------------------- #
#: Radial limits (0-100 scale, pyCirclize convention), innermost to outermost: chord ribbons fill
#: whatever is below IDEOGRAM_R_LIM[0]; the predictor bar track sits in the gap between the two
#: rings; the ideogram (model-colour) ring is outermost, just inside where sector labels land.
#: Cosmetic parameters, tuned by looking at the rendered figure — not a scientific choice.
PREDICTOR_TRACK_R_LIM = (80, 92)
IDEOGRAM_R_LIM = (95, 98)

RIBBON_ALPHA = 0.65

#: Angular gap between adjacent sectors (degrees) — visually separates ribbons leaving neighbouring
#: endpoints, the standard circos look.
SECTOR_SPACE_DEG = 2.0


def _new_polar_axes(cells):
    """A standalone polar-projection figure/axes sized exactly like a normal standalone
    :class:`BasePlot` figure would be (same cells-grid maths ``BasePlot.__init__`` uses) —
    ``stylia.create_figure`` has no ``projection="polar"`` support, and
    ``pycirclize.Circos.plotfig`` requires a real ``PolarAxes``.
    """
    rows, cols = cells
    size = st.get_size()
    width_inches = (cols / CELLS_PER_WIDTH) * size
    height_inches = (rows / CELLS_PER_WIDTH) * size
    fig = plt.figure(figsize=(width_inches, height_inches))
    fig.patch.set_facecolor("white")
    return fig.add_subplot(111, projection="polar")


#: Linkage method for ``order="cluster"`` — UPGMA, the conventional default for this kind of
#: similarity data (no strong reason to prefer complete/ward here; changeable if reviewed).
CLUSTER_LINKAGE_METHOD = "average"


def _ordered_labels(matrix, labels, predictor_label, order, top_n):
    """The ring order for ``labels`` (the bioactivity endpoints, predictor excluded).

    ``order="cluster"``: hierarchical clustering (SciPy's own optimal-leaf-ordering, not a plain
    dendrogram traversal — it minimizes the sum of distances between RING-ADJACENT leaves) on
    ``top_n - shared_count`` as the dissimilarity, so endpoints that share a lot of top-1000
    compounds end up next to each other. Fixes the "looks like an arbitrary mandala" complaint
    structurally: a strong ribbon becomes a short local arc instead of crossing the whole ring, and
    real clusters (or endpoints with nothing to show) become visible at a glance. This generally
    does NOT keep same-model endpoints adjacent, unlike the matrix's own model-grouped order — it
    orders by actual connectivity instead, which is the point.

    ``order`` = anything else (by convention, the pathogen's own ``predictor_family``, e.g.
    ``"abx"``): a direct descending sort by that endpoint's overlap with the predictor node, placed
    clockwise from 12 o'clock (:func:`_node_angles` is already clockwise by construction). Answers a
    different question from clustering — "which endpoints look most like a known antibiotic" — as a
    monotone ranking rather than a structure-revealing layout.
    """
    if order == "cluster":
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        sub = matrix.loc[labels, labels].to_numpy(dtype=float)
        dist = top_n - sub
        np.fill_diagonal(dist, 0.0)
        condensed = squareform(dist, checks=False)
        z = linkage(condensed, method=CLUSTER_LINKAGE_METHOD, optimal_ordering=True)
        return [labels[i] for i in leaves_list(z)]
    return sorted(labels, key=lambda lab: -matrix.loc[lab, predictor_label])


class PathogenChordPlot(BasePlot):
    """Chord (circos-style) diagram, built on pyCirclize: one pathogen's endpoints on a ring, real
    circos ribbons between them (wide at the rim, pinched at the centre — see the module docstring
    for why that needed a library) coloured by their top-1000 shared-actives count, plus a radial
    bar track for their overlap with the pathogen's own predictor node (e.g. ``abx`` for E. coli).

    ``matrix`` and ``model_id`` are the SAME ``sub``/``model_id`` :func:`run_pathogen` already built
    for :class:`PathogenEndpointMatrixPlot` — no new computation happens here. ``predictor_family``
    names the one entry of ``model_id`` that is a predictor rather than a bioactivity endpoint; it is
    excluded from the ribbons and drawn as the bar track instead. ``ax``, if given, must already be a
    ``PolarAxes`` (pyCirclize's own requirement) — the normal call path leaves it ``None`` and lets
    :func:`_new_polar_axes` build one sized to ``cells``.
    """

    def __init__(self, matrix, model_id, predictor_family, ax=None, cells=CELLS,
                 top_n=ACTIVITY_BINARIZE_TOP_N, name=None, order="cluster"):
        if ax is None:
            ax = _new_polar_axes(cells)
        super().__init__(ax=ax, cells=cells)
        self.name = name or "15_endpoint_chord"
        if not matrix.size:
            self._unavailable()
            return

        predictor_idx = model_id.index(predictor_family)
        predictor_label = matrix.columns[predictor_idx]
        label_to_model = dict(zip(matrix.columns, model_id))
        # `natural_labels` fixes the triangularization rule (see below) independent of `order`;
        # `labels` is the actual display order handed to pyCirclize.
        natural_labels = [c for c in matrix.columns if c != predictor_label]
        labels = _ordered_labels(matrix, natural_labels, predictor_label, order, top_n)

        cmap, norm = spectrum_cmap(OVERLAP_MATRIX_SPECTRUM), mcolors.Normalize(vmin=0, vmax=top_n)
        model_colors = _model_color_map(model_id)
        name2hex = {lab: mcolors.to_hex(model_colors[label_to_model[lab]]) for lab in labels}

        # pyCirclize treats its input as a DIRECTED matrix — a full symmetric matrix draws every
        # pair twice (verified against toy data, see module docstring). Triangularized on
        # `natural_labels`' fixed order (not `labels`, which changes with `order`) so the SAME pair
        # is dropped/kept consistently regardless of ring order, then reindexed to display order —
        # reindexing a DataFrame's rows/cols does not change which single cell of a pair is nonzero.
        natural_sub = matrix.loc[natural_labels, natural_labels].astype(float).to_numpy().copy()
        natural_sub[np.tril_indices(len(natural_labels))] = 0.0
        tri = pd.DataFrame(natural_sub, index=natural_labels, columns=natural_labels).loc[labels, labels]

        def link_kws_handler(from_label, to_label):
            v = float(matrix.loc[from_label, to_label])
            return {"fc": mcolors.to_hex(cmap(norm(v))), "alpha": RIBBON_ALPHA, "ec": "none"}

        circos = Circos.chord_diagram(
            tri, space=SECTOR_SPACE_DEG, r_lim=IDEOGRAM_R_LIM, cmap=name2hex, order=labels,
            label_kws=dict(size=st.FONTSIZE_SMALL, color=INK),
            link_kws_handler=link_kws_handler)

        # --- Predictor bar track: one radial bar per endpoint sector, length/colour on the SAME
        # Normalize(0, top_n) scale as the ribbons, so a bar and a ribbon of matching colour
        # represent the same count. Added as a genuine second pyCirclize track per sector, not a
        # hand-drawn patch. A thin INK edge is REQUIRED, not cosmetic: ribbons converge from the
        # ideogram ring right through this track's radius band, and a same-magnitude ribbon is the
        # same colour as the bar sitting on top of it — an unfilled-only ``ec="none"`` bar is
        # verified to disappear into the ribbons behind it (caught by rendering with a forced
        # contrasting edge colour and comparing).
        for lab in labels:
            sector = circos.get_sector(lab)
            track = sector.add_track(PREDICTOR_TRACK_R_LIM)
            v = float(matrix.loc[lab, predictor_label])
            track.bar([sector.size / 2], [v], width=sector.size * 0.85, vmin=0, vmax=top_n,
                      fc=cmap(norm(v)), ec=INK, lw=0.5)

        # pyCirclize defers `colorbar()` (and every other `circos`/track-level convenience call) to
        # a queue executed INSIDE `plotfig()` — called after `plotfig()` has already run, its draw
        # queues but never fires. Verified: moved here (before `plotfig`), it appears; after, it
        # silently does not.
        cb_ax_bounds = (1.08, 0.35, 0.02, 0.35)
        circos.colorbar(bounds=cb_ax_bounds, vmin=0, vmax=top_n, cmap=cmap,
                        label=f"Shared actives (of {top_n})",
                        label_kws=dict(fontsize=st.FONTSIZE_SMALL),
                        tick_kws=dict(labelsize=st.FONTSIZE_SMALL))

        circos.plotfig(ax=self.ax)

        models = list(model_colors)
        self.legend({m: c for m, c in model_colors.items() if m != predictor_family},
                   loc="upper center", bbox_to_anchor=(0.5, -0.05),
                   ncol=min(len(models) - 1, 4))


def save_pathogen_chord_figure(output_dir, matrix, model_id, predictor_family, name,
                               top_n=ACTIVITY_BINARIZE_TOP_N, order="cluster"):
    """Build the chord diagram and merge its footprint into ``figure_cells.json`` (accumulating
    alongside the matrix figure's own entry, same convention as
    :func:`save_pathogen_endpoint_matrix_figure`). ``order`` is passed straight through to
    :class:`PathogenChordPlot` — see :func:`_ordered_labels` for what each mode means.
    """
    plot = PathogenChordPlot(matrix, model_id, predictor_family, name=name, top_n=top_n, order=order)
    footprints = {}
    if plot.is_available:
        plot.save(output_dir)
        footprints[plot.name] = list(plot.cells)
        print(f"  figure: {plot.name}")
    else:
        print(f"  [skip figure] {plot.name}: empty matrix")
    return merge_figure_cells(output_dir, footprints)
