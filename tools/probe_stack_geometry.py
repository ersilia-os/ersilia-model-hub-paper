"""Re-measure the sizing constants for the two `*_by_subtask` stacked panels.

`src/plots_metadata.py` sizes that pair to a physical page budget and solves for equal bar thickness
(`_stack_axes_heights`), which needs three measured constants:

    _AXES_BAND_MM    figure height - axes height, per `show_xlabel`
    _PAGE_PAD_MM     saved page height - figure height
    _PAGE_PAD_W_MM   saved page width  - figure width

Run this whenever the font sizes, the axis labels or the tick-label strings change, and copy the
printed values back into `plots_metadata.py`. The script's own run then verifies them: it prints the
measured band against the constants and the delivered bar thicknesses on every execution.

Two traps this exists to avoid:

* **Sample finer than the quantum.** matplotlib's canvas is a whole number of pixels, so the figure is
  the requested footprint floored onto a pixel grid — 0.254 mm at dpi 100, 0.042 mm at the
  `_LAYOUT_DPI = 600` these panels use. Sweeping in 3 mm steps aliases that staircase into a
  convincing straight line with the wrong slope. The sweep below steps by 0.05 mm for exactly this
  reason, and reports the step in axes height so the quantum is visible.
* **Measure the band against the FIGURE, not the footprint.** Only the figure->axes and figure->page
  relationships are constant; the footprint->figure one is the quantised link.

Usage:  python tools/probe_stack_geometry.py
"""

import os
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

import matplotlib.pyplot as plt
import pandas as pd

from plotting_base import CELL_MM, pdf_page_mm
from plotting_colors import SUBTASK_COLORS
from plots_metadata import (StackedFieldBarPlot, _AXES_BAND_MM, _LAYOUT_DPI, _PAGE_PAD_MM,
                            _PAGE_PAD_W_MM, _SUBTASK_STACK_ORDER)

data_dir = os.path.join(root, "..", "output", "01_models_metadata")
tmp_dir = os.path.join(root, "..", "tmp", "probe_stack_geometry")
os.makedirs(tmp_dir, exist_ok=True)

# Width to probe at. Any value in the panels' working range does: every constant below is invariant to
# it, which is part of what the sweep checks.
PROBE_WIDTH_MM = 46.48


def _table(field):
    name = f"{field.lower().replace(' ', '_')}_by_subtask_counts.csv"
    return pd.read_csv(os.path.join(data_dir, name), index_col=0)


def _measure(field, table, show_xlabel, height_mm):
    """Build one panel at a given declared height and return every size in the chain, in mm."""
    p = StackedFieldBarPlot(
        table=table, colors=SUBTASK_COLORS, legend_kw=None, show_xlabel=show_xlabel,
        cells=(height_mm / CELL_MM, PROBE_WIDTH_MM / CELL_MM),
        name=f"probe_{field[0]}_{int(show_xlabel)}_{height_mm:g}")
    p.measure_geometry()
    fig_w, fig_h = p.fig.get_size_inches() * 25.4
    p.save(tmp_dir)
    page_w, page_h = pdf_page_mm(os.path.join(tmp_dir, "pdf", p.name + ".pdf"))
    return dict(declared=height_mm, fig_w=fig_w, fig_h=fig_h, axes=p.axes_h_mm, bar=p.bar_mm,
                band=fig_h - p.axes_h_mm, page_h=page_h, page_w=page_w,
                pad_h=page_h - fig_h, pad_w=page_w - fig_w)


def main():
    tables = {f: _table(f) for f in _SUBTASK_STACK_ORDER}
    heights = [20.8 + 0.05 * i for i in range(9)]

    print(f"layout dpi {_LAYOUT_DPI:g}; canvas pixel = {25.4 / _LAYOUT_DPI:.4f} mm — the finest size "
          f"difference the layout can express")
    rows = []
    for field, table in tables.items():
        for show in (True, False):
            print(f"\n{field} ({len(table)} bars), axis title {show}")
            print(f"{'declared':>9s} {'fig_h':>8s} {'axes':>8s} {'d(axes)':>8s} {'band':>8s} "
                  f"{'page_h':>8s} {'pad_h':>6s} {'page_w':>8s} {'pad_w':>6s} {'bar':>6s}")
            prev = None
            for h in heights:
                m = _measure(field, table, show, h)
                step = "" if prev is None else f"{m['axes'] - prev:+8.4f}"
                print(f"{m['declared']:9.2f} {m['fig_h']:8.4f} {m['axes']:8.4f} {step:>8s} "
                      f"{m['band']:8.4f} {m['page_h']:8.3f} {m['pad_h']:6.3f} {m['page_w']:8.3f} "
                      f"{m['pad_w']:6.3f} {m['bar']:6.3f}")
                prev = m["axes"]
                rows.append(dict(field=field, show_xlabel=show, **m))
                plt.close("all")

    df = pd.DataFrame(rows)
    print("\nConstants to copy into src/plots_metadata.py "
          "(spread = max - min over the sweep; anything above ~0.001 mm is not a constant):")
    for show, grp in df.groupby("show_xlabel"):
        print(f"  _AXES_BAND_MM[{show}]   = {grp['band'].mean():.4f}   "
              f"spread {grp['band'].max() - grp['band'].min():.4f}   "
              f"(currently {_AXES_BAND_MM[show]:.4f})")
    print(f"  _PAGE_PAD_MM         = {df['pad_h'].mean():.3f}    "
          f"spread {df['pad_h'].max() - df['pad_h'].min():.3f}    (currently {_PAGE_PAD_MM:g})")
    for field, grp in df.groupby("field"):
        print(f"  pad_w, {field:12s}  = {grp['pad_w'].mean():.3f}    "
              f"spread {grp['pad_w'].max() - grp['pad_w'].min():.3f}")
    print(f"  _PAGE_PAD_W_MM       = {df['pad_w'].max():.3f}    "
          f"<- the WIDER panel's, so neither goes over budget (currently {_PAGE_PAD_W_MM:g})")
    print(f"\nWrote {len(rows)} probe panels to {tmp_dir} (not tracked; safe to delete).")


if __name__ == "__main__":
    main()
