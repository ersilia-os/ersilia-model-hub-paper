"""Base plotting module: stylia presets + the BasePlot contract.

This mirrors the plotting architecture of the ``zairachem-docker`` repo
(``zairachem/report/__init__.py``): the publication presets are applied **once at
import time**, and every concrete plot subclasses :class:`BasePlot`.

The key contract is the ``ax`` argument. A plot either

- self-builds a standalone single-panel figure (``ax=None``), or
- draws into an axis handed to it by a figure builder (``ax=<existing axis>``).

The second mode is what lets a builder compose many plots into one condensed
multi-panel figure while each plot class stays reusable on its own.

Semantic colours live in :mod:`plotting_colors`; concrete plots live in the
per-analysis modules (e.g. :mod:`plots_metadata`).
"""

import math
import os

import stylia

from plotting_utils import ref_line as _ref_line
from plotting_utils import swatch_legend as _swatch_legend

# Publication-ready presets, applied once when this module is imported. Any script
# that imports the plotting stack gets the same style without repeating the calls.
stylia.set_format("print")
stylia.set_style("article")

# Reference grid for Nature-style figure footprints. stylia's "print" format is
# SIZE = 7.09" = 180 mm (Nature two-column width), and create_figure's width/height are
# fractions of SIZE. Splitting that width into 6 cells gives true 3 cm square cells, so a
# panel with cells=(rows, cols) renders at (cols/6 x 180 mm) wide by (rows/6 x 180 mm) tall.
# A full page (180 x 215 mm) is therefore 6 cells wide by ~7.2 cells tall.
CELLS_PER_WIDTH = 6
CELL_CM = 3.0


class BasePlot:
    """Base class for a single plot panel.

    Parameters
    ----------
    ax : matplotlib Axes or None
        If ``None``, a standalone one-panel figure is created via
        ``stylia.create_figure``, sized from ``cells`` (see below), and its axis is used.
        Otherwise the given axis is drawn into (used when a builder composes a multi-panel
        figure) and ``cells`` does not affect sizing.
    cells : tuple(int, int)
        Footprint on the reference grid as ``(rows, cols)`` in 3 cm square cells. Drives the
        standalone figure size: ``width = cols / CELLS_PER_WIDTH``, ``height = rows /
        CELLS_PER_WIDTH`` (fractions of stylia's print width).
    """

    def __init__(self, ax=None, cells=(2, 2), **kwargs):
        self.cells = cells
        self.name = "base"
        self.is_available = True
        if ax is None:
            rows, cols = cells
            fig, axs = stylia.create_figure(
                1, 1, width=cols / CELLS_PER_WIDTH, height=rows / CELLS_PER_WIDTH
            )
            self.ax = axs.next()
            self.fig = fig
        else:
            self.ax = ax
            self.fig = ax.figure

    # -- shared chrome helpers (house style flows through these) -------------
    def label(self, xlabel="", ylabel="", title=""):
        """Label this panel's axis via stylia (the one place panels set axis text).

        ``title`` is accepted for call-site compatibility but intentionally NOT rendered: per the
        figure conventions a standalone panel carries no title (titles and lettering are added in
        Illustrator). Small-multiple cell identifiers go through ``stylia.label`` / ``roc_panel``
        directly and are unaffected.
        """
        stylia.label(self.ax, xlabel=xlabel, ylabel=ylabel)

    def ref_line(self, value, axis="x", **kw):
        """Dashed neutral reference line on this panel's axis."""
        return _ref_line(self.ax, value, axis=axis, **kw)

    def legend(self, mapping, **kw):
        """Frameless swatch legend on this panel's axis from a ``{label: colour}`` mapping."""
        return _swatch_legend(self.ax, mapping, **kw)

    def _unavailable(self):
        """Mark the panel as having no data to draw (``save`` becomes a no-op)."""
        self.is_available = False

    def save(self, output_dir):
        """Save this plot as a standalone figure into ``output_dir``.

        Writes both a raster preview (``output_dir/png/<name>.png``, 600 dpi) and a vector
        copy for Illustrator (``output_dir/pdf/<name>.pdf``, fonts embedded), then closes the
        figure. Only meaningful when the plot built its own figure (``ax=None``); a builder
        that composes a shared figure saves the whole figure itself. No-op when the plot is
        not applicable.
        """
        if not self.is_available:
            return
        import matplotlib.pyplot as plt

        png_dir = os.path.join(output_dir, "png")
        pdf_dir = os.path.join(output_dir, "pdf")
        os.makedirs(png_dir, exist_ok=True)
        os.makedirs(pdf_dir, exist_ok=True)
        # Make THIS plot's figure current so stylia.save_figure targets it, regardless of how
        # many other figures are open (plots may be constructed before any are saved).
        plt.figure(self.fig.number)
        stylia.save_figure(os.path.join(png_dir, self.name + ".png"))
        stylia.save_figure(os.path.join(pdf_dir, self.name + ".pdf"))
        plt.close(self.fig)


class MultiPanelPlot(BasePlot):
    """Base for a panel that owns a stylia figure with more than one axis.

    Bypasses ``BasePlot.__init__`` (which makes a 1x1 figure) and builds an
    ``nrows x ncols`` grid sized from the ``cells`` footprint, setting only the attributes
    ``save()`` needs. Subclasses call :meth:`_new_figure` and draw into the returned axes.
    """

    def _new_figure(self, nrows, ncols, cells, name):
        self.cells = cells
        self.name = name
        self.is_available = True
        rows, cols = cells
        fig, axs = stylia.create_figure(
            nrows, ncols, width=cols / CELLS_PER_WIDTH, height=rows / CELLS_PER_WIDTH
        )
        self.fig = fig
        return axs


class GridPlot(BasePlot):
    """Base for a small-multiples grid — one cell per item, drawn by a supplied ``panel_fn``.

    Owns all the grid scaffolding (columns/rows maths, figure creation, edge-only axis
    labels, hiding trailing empty cells) so every grid figure is laid out identically.
    Subclasses prepare ``items`` and call :meth:`build_grid`; ``panel_fn(ax, item, color,
    xlabel, ylabel)`` draws one cell and ``color_fn(item)`` picks its colour.
    """

    def build_grid(self, items, *, cols, name, panel_fn, color_fn=None,
                   edge_xlabel="", edge_ylabel=""):
        n = len(items)
        self.name = name
        self.is_available = n > 0
        if not self.is_available:
            self.cells = (1, cols)
            return
        cols = min(cols, n)
        rows = math.ceil(n / cols)
        self.cells = (rows, cols)
        fig, axs = stylia.create_figure(
            rows, cols, width=cols / CELLS_PER_WIDTH, height=rows / CELLS_PER_WIDTH
        )
        self.fig = fig
        flat = list(axs.axs_flat)
        color_fn = color_fn or (lambda _it: None)
        for i, item in enumerate(items):
            col = i % cols
            xlabel = edge_xlabel if (i + cols) >= n else ""
            ylabel = edge_ylabel if col == 0 else ""
            panel_fn(axs[i], item, color_fn(item), xlabel, ylabel)
        for ax in flat[n:]:
            ax.axis("off")
        self.ax = flat[0]
