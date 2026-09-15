"""Four-track hub timeline on one shared year axis.

Requires data/raw/github_stats/ (run tools/fetch_github_stats.py first) and
output/01_models_metadata/models_over_time_by_task.csv (run 01_ersilia_metadata.py first).

Numbered 01b because it READS 01's output: the Models track is 01's cumulative series, so this
step must run after it. The other three tracks are aggregated here from the GitHub collector's
raw CSVs — this script is the merge of the former steps 08 (community aggregation) and 09 (the
timeline figure), reduced to only what the timeline needs.

EXPLORATORY. Unlike steps 01-07 the GitHub input is not routed through 00_download_data.py — the
collector is a standalone tool because this figure may not survive review. See the header of
tools/fetch_github_stats.py.

Every GitHub count is a snapshot of a moving target. The snapshot date is read from
data/raw/github_stats/snapshot.json, printed, and copied into 01b_snapshot.txt beside the figure
so no caption can be written without it.

The four tracks mix two kinds of quantity, and a caption must say so:
    stocks (models, people)              cumulative — the height is the size of the hub
    flows  (commits, issues)             per month  — the height is a rate
Vertical comparison ACROSS tracks is therefore meaningless. Horizontal comparison is the point.

SCOPE: the three GitHub tracks are scoped to the ~247 repositories that ARE the Model Hub (the
CLI, the packaging / maintenance / workflow infrastructure, the model template, and one repo per
model — see repo_set.csv, written by the collector). That is deliberately NOT the whole
ersilia-os organisation, which also holds websites, grant repos and capstones.

Output
------
    output/01b_community_stats/01b_timeline_series.csv   # every series on one month index
    output/01b_community_stats/01b_snapshot.txt          # GitHub collection date
    output/01b_community_stats/png/hub_timeline.png
    output/01b_community_stats/pdf/hub_timeline.pdf
    output/01b_community_stats/figure_cells.json
"""

import json
import os
import sys

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

from plots_timeline import save_timeline_figure  # noqa: E402
from default import GITHUB_BOT_ACCOUNTS  # noqa: E402

stats_dir = os.path.join(root, "..", "data", "raw", "github_stats")
metadata_dir = os.path.join(root, "..", "output", "01_models_metadata")
outpath = os.path.join(root, "..", "output", "01b_community_stats")
os.makedirs(outpath, exist_ok=True)

with open(os.path.join(stats_dir, "snapshot.json")) as f:
    snapshot = json.load(f)
snapshot_date = snapshot["snapshot_date"]
print(f"GitHub snapshot: {snapshot_date}  (three of four tracks move daily)")

# ---------------------------------------------------------------------------
# Models track — read, not derived
# ---------------------------------------------------------------------------
# The only series this script does not aggregate itself. Step 01 writes it because it is the only
# consumer of Airtable's date columns; reading the summary rather than re-deriving it from the
# frozen snapshot (AIRTABLE_METADATA_FILE in src/default.py) is the repo's summary-CSV rule, and it
# keeps the two scripts from holding two copies of the same date logic.
#
# NOTE its denominator differs from every panel in step 01: the series runs on the UNFILTERED
# metadata (all 224 models on the frozen 2026-08-07 snapshot, whatever their Status), because a model
# in maintenance was still incorporated on its date. Step 01's panels are the 220 Ready models. The
# plotted series reaches 213 rather than 224 because it is trimmed to the last complete month and
# four models carry no Incorporation Date. A caption using this track must quote the plotted number.
models_path = os.path.join(metadata_dir, "models_over_time_by_task.csv")
if not os.path.exists(models_path):
    raise SystemExit(f"Missing {os.path.relpath(models_path, os.path.join(root, '..'))}\n"
                     f"  Run scripts/01_ersilia_metadata.py first.")
models = pd.read_csv(models_path, index_col=0)
models.index = pd.to_datetime(models.index)

# ---------------------------------------------------------------------------
# Hub repo set
# ---------------------------------------------------------------------------
# Membership only — the collector records a core/model group per repo, but every track here is a
# total, so the group is never used. (The former step 08 stacked its commit and issue panels by
# group; those panels are gone.)
repo_set = pd.read_csv(os.path.join(stats_dir, "repo_set.csv"))
hub_repos = set(repo_set["repo"])
print(f"Hub repo set: {len(hub_repos)} repos")

# ---------------------------------------------------------------------------
# People and Commits tracks — from the weekly commit series
# ---------------------------------------------------------------------------
# stats/contributors reports whole WEEKS, so a week is assigned to the month its Monday falls in.
# A handful of commits near a month boundary therefore land one month early. Immaterial at monthly
# resolution, and the alternative — splitting a week's count across two months — would invent
# per-day precision the endpoint does not provide.
#
# Automation accounts are not community members. GitHub reports them all with type "User", so they
# come off an explicit list rather than a field — see GITHUB_BOT_ACCOUNTS.
commits_raw = pd.read_csv(os.path.join(stats_dir, "commit_weeks.csv"))
n_bot_commits = int(commits_raw.loc[commits_raw["login"].isin(GITHUB_BOT_ACCOUNTS), "commits"].sum())
commits_raw = commits_raw[~commits_raw["login"].isin(GITHUB_BOT_ACCOUNTS)]
commits_raw = commits_raw[commits_raw["repo"].isin(hub_repos)]
print(f"Commits: {n_bot_commits:,} bot commits excluded")

# ---------------------------------------------------------------------------
# Issues track — from the org-wide participation sweep, narrowed to the hub
# ---------------------------------------------------------------------------
# org_participation.csv records both issues and pull requests across the whole organisation. Two
# narrowings, both deliberate: to the hub repo set (the sweep is org-wide), and to issues only —
# a PR is a different act, and the track is labelled "Issues".
part = pd.read_csv(os.path.join(stats_dir, "org_participation.csv"))
part = part[~part["login"].isin(GITHUB_BOT_ACCOUNTS)]
part["repo_name"] = part["repo"].str.split("/").str[-1]
hub_part = part[part["repo_name"].isin(hub_repos)].copy()
hub_issues = hub_part[hub_part["kind"] == "issue"]


def monthly_total(frame, date_col, start=None, end=None):
    """Rows per month, reindexed onto a gapless month range ending at ``month_end`` timestamps.

    ``start``/``end`` force a shared range across series that are read against each other —
    without them, commits (from 2020) and issues (from 2022) would sit on different time bases and
    their features would not line up under the shared x axis.
    """
    month = pd.to_datetime(frame[date_col]).dt.to_period("M")
    idx = pd.period_range(start or month.min(), end or month.max(), freq="M")
    series = month.value_counts().reindex(idx, fill_value=0).sort_index()
    series.index = idx.to_timestamp(how="end")
    return series


def cumulative_first_seen(frame, date_col):
    """Cumulative distinct logins, each counted from the month of their FIRST appearance.

    Counting first appearances rather than activity is what makes this a count of people: someone
    who comes back every month is one person, not twelve.
    """
    month = pd.to_datetime(frame[date_col]).dt.to_period("M")
    idx = pd.period_range(month.min(), month.max(), freq="M")
    first = month.groupby(frame["login"]).min()
    series = first.value_counts().reindex(idx, fill_value=0).sort_index().cumsum()
    series.index = idx.to_timestamp(how="end")
    return series


# A weekly row carries a commit COUNT, so it is expanded to one row per commit before binning —
# otherwise the monthly figure would count weeks, not commits.
commits_expanded = commits_raw.loc[commits_raw.index.repeat(commits_raw["commits"])]

# Commits and issues share one time base: the first month either series has data onwards.
span_start = min(pd.to_datetime(commits_expanded["week"]).dt.to_period("M").min(),
                 pd.to_datetime(hub_issues["created_at"]).dt.to_period("M").min())
commits = monthly_total(commits_expanded, "week", span_start)
issues = monthly_total(hub_issues, "created_at", span_start)

# Two ways of being a contributor to the hub's own code. Only "Commit authors" is PLOTTED (see the
# tracks list below); the issue/PR series is carried into the CSV because the gap between them is
# the caveat any caption about "contributors" needs — 107 people wrote code, 335 took part.
commit_authors = cumulative_first_seen(commits_raw, "week")
issue_pr_authors = cumulative_first_seen(hub_part, "created_at")

# ---------------------------------------------------------------------------
# Trim the incomplete snapshot month, for plotting only
# ---------------------------------------------------------------------------
# The month the snapshot was taken is INCOMPLETE. Plotted as though it were a whole month it is a
# drop to zero at the right edge of three tracks, which reads as the project stopping. The figure
# therefore stops at the last COMPLETE month. Nothing is discarded: the CSV below keeps the
# partial month, and its counts are printed.
last_complete = pd.Period(snapshot_date, freq="M") - 1


def trim(series):
    return series[series.index.to_period("M") <= last_complete]


n_partial_commits = int(commits[commits.index.to_period("M") > last_complete].sum())
n_partial_issues = int(issues[issues.index.to_period("M") > last_complete].sum())
print(f"Partial month {snapshot_date[:7]} excluded from the plotted series "
      f"(kept in the CSV): {n_partial_commits} commits, {n_partial_issues} issues")

models_total = trim(models.sum(axis=1))
commits_plot = trim(commits)
issues_plot = trim(issues)
# Cumulative integer counts of people: a leading gap is "nobody yet" (0) and an interior gap
# carries the last value forward. Never interpolated — a fractional person is not a thing, and a
# straight line between two months would invent arrivals.
people_plot = trim(commit_authors).fillna(0).cummax()

# One shared x range across all four tracks, so every track is drawn on the same axis even though
# models start 2020-11 and the first commit lands 2020-06.
_all = (models_total, people_plot, commits_plot, issues_plot)
xlim = (min(s.index.min() for s in _all), max(s.index.max() for s in _all))
print(f"Shared x axis: {xlim[0].date()} to {xlim[1].date()}")

# Workflow runs are deliberately NOT a track. GitHub deletes run records after ~13-14 months, so
# that series only covers 2025-06 onward; on a six-year shared axis it is empty for three quarters
# of its width, which costs a whole track to say almost nothing. The collector no longer fetches
# them.
tracks = [
    {"label": "Models", "series": models_total, "kind": "stock"},
    {"label": "People", "series": people_plot, "kind": "stock"},
    {"label": "Commits", "series": commits_plot, "kind": "flow"},
    {"label": "Issues", "series": issues_plot, "kind": "flow"},
]

# One CSV carrying every series on one month index, so the whole figure is reproducible from a
# single file. UNTRIMMED — the partial snapshot month is kept here even though it is not plotted.
# This file is now the only on-disk record of these series (the former step 08 wrote four CSVs that
# each kept the partial month), so dropping it here would lose it entirely.
#
# The series start in different months, so the joined frame has NaNs at the head of some columns.
# For the two people columns those are filled: they are cumulative counts of people, so a month
# before the first arrival is a real 0 ("nobody yet"), not missing data, and any interior gap
# carries the last value forward. Never interpolated — a fractional person is not a thing, and a
# straight line between two months would invent arrivals.
#
# A TRAILING NaN is left as NaN and means something different: that source series ends before the
# snapshot month. `models_cumulative` stops at the last month carrying an Incorporation Date, so it
# is blank for the partial month — which is not the same claim as "no models were added".
people = pd.DataFrame({"commit_authors_cumulative": commit_authors,
                       "issue_pr_authors_cumulative": issue_pr_authors})
people = people.apply(lambda s: s.fillna(0).cummax() if s.first_valid_index() is not None else s)
combined = pd.DataFrame({
    "models_cumulative": models.sum(axis=1),
    "commit_authors_cumulative": people["commit_authors_cumulative"],
    "issue_pr_authors_cumulative": people["issue_pr_authors_cumulative"],
    "commits_per_month": commits,
    "issues_per_month": issues,
})
combined.index.name = "month_end"
combined.to_csv(os.path.join(outpath, "01b_timeline_series.csv"))

# Two sets of totals, because they differ and a caption must use the right one: the ALL-TIME
# figures include the incomplete snapshot month, the PLOTTED figures stop at the last complete
# month. At the 2026-08-02 snapshot that is 336 vs 335 issue/PR authors and 1,470 vs 1,469 issues —
# small, but a reader checking the figure against this file would otherwise find a mismatch and no
# explanation for it.
with open(os.path.join(outpath, "01b_snapshot.txt"), "w") as f:
    f.write(f"GitHub snapshot date: {snapshot_date}\n"
            f"Organisation: {snapshot['org']} "
            f"({snapshot['org_meta']['public_repos']} public repos)\n"
            f"Hub repo set: {len(hub_repos)} repositories\n"
            f"\n"
            f"ALL TIME (includes the incomplete snapshot month {snapshot_date[:7]}):\n"
            f"  Distinct hub commit authors:      {int(commit_authors.iloc[-1])}\n"
            f"  Distinct hub issue / PR authors:  {int(issue_pr_authors.iloc[-1])}\n"
            f"  Commits:                          {int(commits.sum())}\n"
            f"  Issues:                           {int(issues.sum())}\n"
            f"\n"
            f"AS PLOTTED (through the last complete month, {last_complete}):\n"
            f"  Models (cumulative):              {int(models_total.iloc[-1])}\n"
            f"  People / commit authors:          {int(people_plot.iloc[-1])}\n"
            f"  Commits:                          {int(commits_plot.sum())}\n"
            f"  Issues:                           {int(issues_plot.sum())}\n"
            f"\n"
            f"The Models track runs on the UNFILTERED metadata and is trimmed to the last complete "
            f"month, so it reaches {int(models_total.iloc[-1])} where step 01's panels count only "
            f"Ready models on the frozen snapshot. The two numbers are not meant to agree.\n"
            f"The three GitHub tracks move daily; a caption is only reproducible against the "
            f"snapshot date above.\n"
            f'The People track counts COMMIT AUTHORS only. "Contributors" in the wider sense '
            f"(issue / PR authors) is roughly three times larger — see the two counts above.\n")

save_timeline_figure(tracks=tracks, xlim=xlim, output_dir=outpath)
