"""The hub timeline: four thin tracks on one shared year axis (script 01b).

This is the repo's only time-series figure. The four series it draws — models, people, commits,
issues — are aggregated by ``scripts/01b_community_stats.py``, except the Models track, which it
reads from ``output/01_models_metadata/models_over_time_by_task.csv``. Drawing them as four
separate panels was the alternative and does not work: each would carry its own x axis, so nothing
could be read across them, and they cannot simply be stacked afterwards because
``stylia.save_figure`` crops every file to its own content and their left edges land at different
x positions.

This module draws every series into **one** figure with stacked axes, which is what
``docs/figure_conventions.md`` already prescribes for exactly this problem: describing the 0.25 mm
mis-registration between the two stacked subtask panels, it says "fixing it properly means drawing
both into one figure with a shared x axis (``MultiPanelPlot``)". So this is not an exception to
the one-chart-per-file rule; it is the documented remedy for a shared axis.

Two kinds of quantity share the figure and a caption must say so:

- **stocks** (models, people) are cumulative — the height *is* the size of the hub;
- **flows** (commits, issues) are per-month counts — the height is a rate.

They are not comparable to each other vertically. What IS comparable is the horizontal position
of every feature, which is the entire point of the shared axis.

The tracks are only ~4.4 mm tall, and nearly every layout choice here follows from that: axes
placed by hand rather than as subplots, two y ticks each, horizontal track labels, and vertical
gridlines only. Each is explained where it is set.
"""

import warnings

import matplotlib.dates as mdates
import matplotlib.ticker as mticker
import numpy as np
import stylia
from matplotlib.lines import Line2D

from plotting_base import MultiPanelPlot
from plotting_colors import (REFERENCE_LINE, TIMELINE_COLORS, TIMELINE_FILL_ALPHA,
                             TIMELINE_HUES, TIMELINE_SECONDARY_LIGHTEN, hue)
from plotting_utils import LEGEND_KW

#: Full page width, four very thin tracks: 29.4 mm nominal, which crops to **32.3 mm** —
#: the size actually placed on the page, and the number that was tuned. The nominal figure is
#: smaller than the crop because the outermost tick labels sit outside the axes rectangle; picking
#: `rows` from the nominal height would miss the target by ~10%.
#:
#: Each track gets ~4.4 mm of drawing height. Everything below (the inter-track gap, which ticks get
#: labels, hand-placed axes) is a consequence of that budget. The script prints the measured
#: per-track height on every run so the margin stays checkable rather than assumed.
TIMELINE_CELLS = (0.98, 6)

#: One year tick per year. The 45 mm community panels need ``YearLocator(2)`` because six labels
#: collide at that width; at 180 mm all seven fit with room to spare.
_YEAR_TICK_STEP = 1

#: Gap between tracks, in figure coordinates. Not cosmetic: with the tracks touching, a track's "0"
#: label and the label at the top of the track below it are centred on the same spine and overlap.
#: The floor is therefore one tick-label height; this is that, measured (see ``_place_axes``), plus
#: a little air. Every mm here comes straight out of the tracks, so it is kept as small as the
#: labels allow rather than set by eye.
_TRACK_GAP = 0.06

#: Axes rectangle in figure coordinates. The tracks are placed by hand inside it (see
#: ``_place_axes``).
#:
#: Left reserves the horizontal track label plus the widest tick label ("1,000"); bottom reserves
#: the single shared year axis. Both are deliberately a little generous: ``bbox_inches="tight"``
#: crops whatever is unused, so over-reserving costs nothing on the page, while under-reserving
#: clips a label.
_MARGIN_LEFT = 0.105
_MARGIN_RIGHT = 0.995
_MARGIN_BOTTOM = 0.20
_MARGIN_TOP = 0.98

#: Bins handed to ``MaxNLocator`` purely to pick a ROUND top value — 160 rather than 215, 1,000
#: rather than 1,366. Only the first and last of the positions it returns survive
#: (see ``_thin_y_labels``), so this does not control how many ticks are drawn; two always are.
_Y_TICK_BINS = 3

#: x position of the horizontal track label, in axes coordinates, i.e. this far LEFT of the axes.
#: Set past the widest tick label ("1,000") so the two never collide. Because every track uses the
#: same value, the four labels are left-aligned with each other by construction — which is also
#: what ``fig.align_ylabels`` would do, kept as well since it costs nothing.
_YLABEL_X = -0.055


def _line_handle(color, label):
    """Legend handle for one track line, styled like the line actually plotted."""
    return Line2D([], [], color=color, linewidth=1.0, label=label)


class HubTimelinePlot(MultiPanelPlot):
    """Stacked thin tracks sharing one x axis.

    ``tracks`` is an ordered list of dicts, each with:

    ``label``   y-axis label (also the ``TIMELINE_COLORS`` key)
    ``series``  a Series (single line) or dict of ``{name: Series}`` (multiple lines), month index
    ``kind``    ``"stock"`` (cumulative) or ``"flow"`` (per-month) — recorded for the caption and
                used to pick the y tick strategy, not to transform the data
    ``no_data_before``  optional Timestamp; everything left of it is shaded as unavailable

    ``xlim`` is applied identically to every axis. ``MultiPanelPlot`` does not set ``self.ax``,
    so the ``BasePlot`` chrome helpers are unavailable here and labelling goes through
    ``stylia.label`` directly — the same way the five subclasses in ``plots_chembl_curation``
    do it.
    """

    def __init__(self, tracks, xlim, cells=TIMELINE_CELLS, name="hub_timeline"):
        # _new_figure builds an nrows x 1 grid of SUBPLOTS, and subplots are exactly what must be
        # avoided here — see _place_axes. Ask it for a 1 x 1 figure purely to get the right canvas
        # size and the bookkeeping (name, cells, fig), then throw the axis away and place the
        # tracks by hand.
        self._new_figure(1, 1, cells, name)
        self.is_available = bool(tracks)
        if not self.is_available:
            return
        for stray in list(self.fig.axes):
            stray.remove()

        self.tracks = tracks
        self.axes = self._place_axes(len(tracks))
        self.totals = {}

        for i, (ax, track) in enumerate(zip(self.axes, tracks)):
            self._track(ax, track, xlim, is_last=(i == len(tracks) - 1))

        # Without this the y labels sit at different x positions, because each is placed relative
        # to its own tick-label column and "1,000" is wider than "200" — which made the Commits
        # label visibly stick out to the left of the other three. align_ylabels pins them all to
        # the leftmost requirement.
        self.fig.align_ylabels(self.axes)
        self._thin_y_labels()

    def _place_axes(self, n):
        """Lay the tracks out by hand, with an exact gap between them and no gridspec.

        The gap between tracks has to be set exactly, and with subplots it cannot be set at all.
        ``stylia.save_figure`` calls ``plt.tight_layout()`` unconditionally, tight_layout separates
        gridspec rows by ``h_pad`` (default 1.08 font-size units, ~2.3 mm here), and it recomputes
        from the gridspec every time — so ``subplots_adjust(hspace=...)`` is silently discarded and
        there is no argument to pass through stylia's wrapper.

        ``add_axes`` creates axes with no subplotspec, which tight_layout leaves alone. That gives
        exact control of the gap, and returns the padding tight_layout was spending to the tracks.
        """
        left, right = _MARGIN_LEFT, _MARGIN_RIGHT
        bottom, top = _MARGIN_BOTTOM, _MARGIN_TOP
        height = (top - bottom - _TRACK_GAP * (n - 1)) / n
        # Built top-down so axes[0] is the top track, matching the order tracks are passed in.
        return [self.fig.add_axes([left, top - (i + 1) * height - i * _TRACK_GAP,
                                   right - left, height])
                for i in range(n)]

    def _thin_y_labels(self):
        """Label only the bottom and top gridline of each track.

        One draw for the whole figure, not one per track: the locators have to be realised before
        their positions can be read back, and doing that inside the track loop meant four full
        draws and four copies of every layout warning.

        Only ticks inside the view can be labelled. ``MaxNLocator`` routinely places one above the
        data — 0/80/160/240 for a 215 maximum — so taking the raw last entry labelled a tick that
        is never drawn, which is why the top label went missing on the first attempt.
        """
        self.fig.canvas.draw()
        for ax in self.axes:
            lo, hi = ax.get_ylim()
            visible = [t for t in ax.get_yticks() if lo <= t <= hi]
            # Ticks are REMOVED, not blanked. While horizontal gridlines were drawn, the
            # intermediate ticks earned their place by carrying a rule across the track; with the
            # rules gone they would leave bare unlabelled marks on the spine. The locator still
            # runs first, so the top tick is a round number chosen by MaxNLocator rather than the
            # raw data maximum.
            ax.set_yticks([visible[0], visible[-1]] if len(visible) > 1 else visible)
            ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _pos: f"{v:,.0f}"))

    # -- one track ----------------------------------------------------------
    def _track(self, ax, track, xlim, is_last):
        label = track["label"]
        base = TIMELINE_COLORS[label]
        series = track["series"]
        if not isinstance(series, dict):
            series = {label: series}

        # More than one line in a track means one quantity measured two ways, so the extra series
        # take lighter weights of the SAME hue — a second hue would read as a second category and
        # break the one-hue-per-track scheme.
        #
        # The weights go by SIZE, not by the order the caller listed them, and the series are
        # drawn largest first. These curves nest (every commit author is also an issue/PR author),
        # so the small one sits inside the big one's filled area: give the big one the base hue and
        # the small one the tint and the small one vanishes into the wash. Largest therefore takes
        # the palest weight and is drawn at the back; smallest takes the full hue on top.
        order = sorted(series, key=lambda n: float(np.nanmax(series[n].values)), reverse=True)
        tint = hue(TIMELINE_HUES[label], lighten=TIMELINE_SECONDARY_LIGHTEN)
        weights = {n: (tint if i < len(order) - 1 else base) for i, n in enumerate(order)}

        for name in order:
            s = series[name]
            color = weights[name]
            x = s.index.to_pydatetime()
            y = np.asarray(s.values, dtype=float)
            ax.fill_between(x, y, color=color, alpha=TIMELINE_FILL_ALPHA, linewidth=0)
            ax.plot(x, y, color=color, linewidth=1.0)
            self.totals[name] = {"last": float(y[-1]) if len(y) else float("nan"),
                                 "sum": float(np.nansum(y)),
                                 "kind": track.get("kind", "stock")}

        # Region where the data does not exist, as opposed to being zero. Drawn UNDER the series
        # and across the full track height, with the note inside the band: an empty stretch of
        # track is otherwise indistinguishable from a stretch of no activity.
        edge = track.get("no_data_before")
        if edge is not None:
            ax.axvspan(xlim[0], edge, facecolor=REFERENCE_LINE, alpha=0.25,
                       linewidth=0, zorder=0)
            ax.text(xlim[0] + (edge - xlim[0]) / 2, 0.5, "no data retained",
                    transform=ax.get_xaxis_transform(), ha="center", va="center",
                    fontsize=stylia.FONTSIZE_SMALL, color=REFERENCE_LINE, zorder=1)

        ax.set_xlim(*xlim)
        ax.set_ylim(0, None)
        # Gridlines and tick LABELS are separate requirements here and are satisfied separately.
        # The locator puts four rules across the track so the eye has a scale to measure spikes
        # against; the formatter then blanks everything but the floor and the top, because four
        # labels at ~2.1 mm of type each need ~8 mm against a ~5 mm track and would overlap into
        # mush. The rules give the scale, the two end labels give the range.
        ax.yaxis.set_major_locator(mticker.MaxNLocator(_Y_TICK_BINS, integer=True))
        ax.xaxis.set_major_locator(mdates.YearLocator(_YEAR_TICK_STEP))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y"))
        # Vertical rules only. The year lines are what the shared axis is for — they run down the
        # whole figure and are how a feature in one track is located against another. Horizontal
        # rules cross them every few mm in a track this short, and the resulting mesh reads as
        # texture rather than as a scale; each track's two y labels carry the range instead.
        ax.yaxis.grid(False)
        ax.xaxis.grid(True)
        if not is_last:
            # Marks AND labels go on every track but the bottom one. The gap between tracks is
            # only ~1.8 mm, so an upper track's x ticks would hang most of the way into the track
            # below and read as marks on that track's data rather than as its own axis. The year
            # gridlines already carry the registration the ticks used to provide, and one shared
            # axis at the foot of the figure is the whole premise. (Stricter than
            # StackedFieldBarPlot.show_x, which keeps its marks because its panels sit further
            # apart.)
            ax.tick_params(axis="x", bottom=False, labelbottom=False)
        stylia.label(ax, xlabel="", ylabel=label)
        # HORIZONTAL track label, against the house default of a rotated one. Rotated, "Commits"
        # is ~7 mm of type standing in a ~4 mm track: the four labels overlapped each other into a
        # single unreadable column ("IssuesCommitsPeopleModels"). Laid flat, each is ~2.1 mm tall
        # and fits, at the cost of ~8 mm of the 180 mm width — 4%, which this figure can afford
        # and the vertical budget cannot.
        ax.yaxis.label.set(rotation=0, ha="right", va="center")
        ax.yaxis.set_label_coords(_YLABEL_X, 0.5, transform=ax.transAxes)

        if len(series) > 1:
            # Legend in the caller's order, which is the order that reads as a sentence, not the
            # size order the drawing used.
            ax.legend(handles=[_line_handle(weights[n], n) for n in series],
                      loc="upper left", **LEGEND_KW)

    # -- reporting ----------------------------------------------------------
    def track_height_mm(self):
        """Measured drawing height of one track, for the legibility check the script prints."""
        self.fig.canvas.draw()
        box = self.axes[0].get_window_extent()
        return box.height / self.fig.dpi * 25.4


def save_timeline_figure(tracks, xlim, output_dir, cells=TIMELINE_CELLS):
    """Render the timeline and write ``figure_cells.json`` for ``output_dir``.

    Sole writer of its manifest, so it writes wholesale rather than merging into an existing one.
    """
    import json
    import os

    # A regression guard, not a workaround for a live problem. At 31.5 mm this figure sits close to
    # the height where matplotlib's tight_layout gives up ("cannot make Axes height small enough to
    # accommodate all Axes decorations") — with every gridline labelled it DID fail, and only
    # blanking the two intermediate labels in _thin_y_labels bought back enough room. Since
    # ``stylia.save_figure`` calls ``tight_layout()`` unconditionally, the warning is caught here
    # and reported as one plain line, so that shrinking the figure or adding a decoration says so
    # instead of scrolling a matplotlib warning past. If it ever fires: the axes stay in one
    # gridspec column and so remain aligned, and bbox_inches="tight" still crops correctly, but the
    # margins revert to matplotlib defaults rather than fitted ones.
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        plot = HubTimelinePlot(tracks, xlim=xlim, cells=cells)
        height_mm = plot.track_height_mm()
        plot.save(output_dir)
    messages = [str(w.message) for w in caught]
    # "not compatible with tight_layout" is EXPECTED and is the mechanism working: the tracks are
    # placed with add_axes precisely so tight_layout leaves them alone. Swallowed silently.
    # "cannot make Axes height small enough" is the real signal — it means the figure has been
    # shrunk or decorated past what its labels need — so only that one is surfaced.
    overconstrained = any("make Axes height small enough" in m for m in messages)

    with open(os.path.join(output_dir, "figure_cells.json"), "w") as f:
        json.dump({plot.name: list(plot.cells)}, f, indent=2)

    print(f"\n[{plot.name}] {len(tracks)} tracks over "
          f"{xlim[0].date()} to {xlim[1].date()}, one shared x axis")
    pad = " " * (len(plot.name) + 3)
    print(f"{pad}measured track height {height_mm:.1f} mm; two y ticks per track, "
          f"vertical year gridlines only")
    if overconstrained:
        print(f"{pad}WARNING tight_layout gave up — the figure is now too short for its "
              f"decorations. Axes stay aligned and the crop is still correct, but margins are "
              f"matplotlib defaults rather than fitted. Raise the height or drop a decoration.")
    # A stock's headline number is where it ends; a flow's is its total, because the last month of
    # a rate series is just one month and says nothing about the whole.
    for name, t in plot.totals.items():
        summary = (f"reaches {t['last']:,.0f}" if t["kind"] == "stock"
                   else f"{t['sum']:,.0f} in total, last plotted month {t['last']:,.0f}")
        print(f"{'':{len(plot.name) + 3}}{name}: {summary}")
    return plot
