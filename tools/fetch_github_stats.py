"""Collect the GitHub series behind the hub timeline (``scripts/01b_community_stats.py``).

DEVIATION FROM THE REPO CONVENTION, ON PURPOSE. Every other external data dependency is
routed through ``scripts/00_download_data.py``. This one is standalone because the timeline
figure is exploratory and may be dropped. If it is kept, fold this into a "Section 5 — GitHub"
of that script and delete this file; the three writers below are already shaped like the
declarative sources in Section 1.

Scoped to exactly what the timeline needs: issue/PR authors, the hub repo set, and weekly
commits. Three further collectors (self-reported contributor locations, per-contributor weeks on
the main repo, and workflow runs) were removed on 2026-08-04 with the community panels they fed —
see git history if any of them is wanted back. Files they already wrote are left in
``data/raw/github_stats/`` (raw data is never cleaned) but are no longer refreshed or recorded in
``snapshot.json``.

Writes only small **summary** CSVs into ``data/raw/github_stats/`` — never per-commit dumps.
Each run stamps ``snapshot.json`` with the date and the headline counts, because unlike the
Airtable export these numbers move every day and a figure built from them is meaningless
without the date it was taken.

Requires an authenticated ``gh`` CLI (``gh auth status``). Read-only throughout.

Usage
-----
    python tools/fetch_github_stats.py           # skip files that already exist
    python tools/fetch_github_stats.py --refresh # re-fetch everything
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from datetime import date

import pandas as pd

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))

ORG = "ersilia-os"
MAIN_REPO = f"{ORG}/ersilia"

outdir = os.path.join(root, "..", "data", "raw", "github_stats")
os.makedirs(outdir, exist_ok=True)

# ---------------------------------------------------------------------------
# The Ersilia Model Hub repo set
# ---------------------------------------------------------------------------
# The hub is not one repository. It is the CLI plus the packaging, maintenance and workflow
# infrastructure around it, plus one repository per model. Anything scoped to ersilia-os/ersilia
# alone measures a fraction of the work; anything scoped to the whole org sweeps in unrelated
# projects (websites, grant repos, capstones). This is the middle scope: the software that makes
# the hub run, and the models it serves.
CORE_REPOS = [
    "ersilia",                  # the CLI and hub itself
    "ersilia-pack",             # model packaging
    "ersilia-pack-utils",       # packaging helpers
    "ersilia-maintenance",      # scheduled model health checks
    "ersilia-model-workflows",  # CI workflows shared by the model repos
    "eos-template",             # the template every model repo is created from
    "ersilia-apptainer",        # Apptainer/Singularity images
]

# Model repositories are named ``eos`` + exactly four alphanumerics. The prefix alone is not
# enough: eos-template, eosvc, eosbench, eosframes, eosquality, eosdev and several others share it
# without being models. eos-template is in CORE_REPOS above precisely because it is not a model.
MODEL_REPO_RE = r"^eos[0-9a-z]{4}$"

# The org's first repo predates the model hub; the search sweep starts here.
FIRST_YEAR = 2020

# GitHub's search API caps any single query at 1,000 results however you paginate, and the org
# has ~3,300 issues. Queries are therefore sliced by quarter, which keeps every slice two orders
# of magnitude below the cap. Slicing is a pagination workaround, not a filter — every issue in
# the range is still collected.
QUARTERS = [("01-01", "03-31"), ("04-01", "06-30"),
            ("07-01", "09-30"), ("10-01", "12-31")]

# The search endpoint is rate-limited to 30 requests/minute for authenticated users, well below
# the 5,000/hour core limit. A short sleep between search calls keeps the sweep under it without
# needing retry-after handling.
SEARCH_PAUSE_S = 2.5

# stats/* endpoints return HTTP 202 with an empty body while GitHub computes the series, then the
# real payload on a later call. Observed live on this repo, so it is handled rather than hoped for.
STATS_RETRIES = 5
STATS_PAUSE_S = 3.0


def gh(*args):
    """Run a read-only ``gh`` command and return its parsed JSON stdout."""
    proc = subprocess.run(["gh", *args], capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"gh {' '.join(args)} failed:\n{proc.stderr.strip()}")
    return json.loads(proc.stdout or "null")


def _write(name, frame):
    path = os.path.join(outdir, name)
    frame.to_csv(path, index=False)
    print(f"    wrote {name}  ({len(frame)} rows)")
    return path


def _skip(name, refresh):
    path = os.path.join(outdir, name)
    if os.path.exists(path) and not refresh:
        print(f"    {name} exists, skipping (use --refresh to re-fetch)")
        return True
    return False


# ---------------------------------------------------------------------------
# 1. Org-wide issue and PR authors
# ---------------------------------------------------------------------------
def fetch_org_participation(refresh):
    """One row per issue/PR across the whole org: author, kind, repo, creation date.

    Commits alone undercount a community badly — most first contributions are an issue or a
    review comment, and the model repos each have their own contributor set that never touches
    the main repo. This is the only sweep that sees all 421 repos.
    """
    name = "org_participation.csv"
    if _skip(name, refresh):
        return
    print("  org-wide issues and PRs (quarterly sweep)...")
    rows = []
    for year in range(FIRST_YEAR, date.today().year + 1):
        for start, end in QUARTERS:
            window = f"{year}-{start}..{year}-{end}"
            for kind in ("issue", "pr"):
                # --slurp wraps the paginated pages into one array. It cannot be combined with
                # --jq (gh rejects the pair), so the fields are picked out here instead.
                pages = gh("api", "-X", "GET", "search/issues",
                           "-f", f"q=org:{ORG} is:{kind} created:{window}",
                           "-f", "per_page=100", "--paginate", "--slurp")
                time.sleep(SEARCH_PAUSE_S)
                for page in (pages or []):
                    for it in page.get("items", []):
                        user = it.get("user") or {}
                        if not user.get("login"):     # deleted account
                            continue
                        rows.append({
                            "login": user["login"], "kind": kind,
                            "created_at": it["created_at"][:10],
                            # repo is not returned directly; recover it from the HTML URL
                            "repo": "/".join(it["html_url"].split("/")[3:5]),
                        })
            print(f"    {window}: {len(rows)} cumulative")
    _write(name, pd.DataFrame(rows).drop_duplicates().sort_values("created_at"))


# ---------------------------------------------------------------------------
# 2. The repo set: core infrastructure + one repo per model
# ---------------------------------------------------------------------------
def fetch_repo_set(refresh):
    """Resolve CORE_REPOS + every model repo to ``repo, group, created_at, archived``."""
    name = "repo_set.csv"
    if _skip(name, refresh):
        return pd.read_csv(os.path.join(outdir, name))
    print("  resolving repo set...")
    listing = gh("repo", "list", ORG, "--limit", "1000",
                 "--json", "name,createdAt,isArchived")
    by_name = {r["name"]: r for r in listing}
    rows = []
    for repo in CORE_REPOS:
        if repo not in by_name:
            raise RuntimeError(f"{ORG}/{repo} not found — check CORE_REPOS")
        rows.append({"repo": repo, "group": "core"})
    for repo in by_name:
        if re.fullmatch(MODEL_REPO_RE, repo):
            rows.append({"repo": repo, "group": "model"})
    for row in rows:
        meta = by_name[row["repo"]]
        row["created_at"] = meta["createdAt"][:10]
        row["archived"] = bool(meta["isArchived"])
    frame = pd.DataFrame(rows).sort_values(["group", "repo"])
    n_model = int((frame["group"] == "model").sum())
    print(f"    {len(frame)} repos: {len(CORE_REPOS)} core + {n_model} model "
          f"({int(frame['archived'].sum())} archived)")
    _write(name, frame)
    return frame


# ---------------------------------------------------------------------------
# 3. Weekly commits per author, per repo, across the whole set
# ---------------------------------------------------------------------------
def fetch_commit_weeks(repos, refresh):
    """``repo, login, week, commits`` for every repo in the set.

    ``stats/contributors`` answers three questions in one call per repo — who has committed, how
    many commits, and in which week — so this is the source for both the commit series and the
    ecosystem contributor count.

    The endpoint returns HTTP 202 with an empty body the first time it is asked about a repo,
    while GitHub computes the series in the background. With ~250 repos, retrying each one five
    times in place would waste most of the run waiting. Instead every repo is asked once per
    round, whatever is ready is banked, and only the stragglers go into the next round — the
    first round doubles as the "warm the cache" pass.

    A repo that is still empty after the last round is recorded with zero commits rather than
    dropped, because an empty result is also what a genuinely commit-less repo returns and the
    two cannot be told apart from the response alone.
    """
    name = "commit_weeks.csv"
    if _skip(name, refresh):
        return
    print(f"  weekly commits for {len(repos)} repos...")
    rows, pending = [], list(repos)
    for rnd in range(STATS_RETRIES):
        still = []
        for repo in pending:
            data = gh("api", f"repos/{ORG}/{repo}/stats/contributors")
            if not data:
                still.append(repo)
                continue
            for entry in data:
                login = (entry.get("author") or {}).get("login")
                if login is None:            # deleted account; commits survive, identity does not
                    continue
                for w in entry["weeks"]:
                    if w["c"]:
                        rows.append({"repo": repo, "login": login,
                                     "week": w["w"], "commits": w["c"]})
        print(f"    round {rnd + 1}: {len(pending) - len(still)} resolved, "
              f"{len(still)} still computing")
        pending = still
        if not pending:
            break
        time.sleep(STATS_PAUSE_S * 2)
    if pending:
        print(f"    {len(pending)} repos returned no data after {STATS_RETRIES} rounds "
              f"(empty repo or stats still cold): {pending[:10]}")
    frame = pd.DataFrame(rows)
    frame["week"] = pd.to_datetime(frame["week"], unit="s").dt.date
    print(f"    {frame['login'].nunique()} distinct commit authors, "
          f"{int(frame['commits'].sum())} commits, "
          f"{frame['repo'].nunique()}/{len(repos)} repos with commits")
    _write(name, frame.sort_values(["week", "repo", "login"]))


# ---------------------------------------------------------------------------
# 4. Snapshot record
# ---------------------------------------------------------------------------
def write_snapshot():
    """Date-stamp the collection, with the headline counts these CSVs imply.

    Every number here moves daily. A figure built from this directory is only reproducible
    against the date recorded in this file, which is why it is written on every run even when
    the CSVs themselves were skipped.
    """
    repo = gh("api", f"repos/{MAIN_REPO}",
              "--jq", "{stars: .stargazers_count, forks: .forks_count, "
                      "created_at, pushed_at}")
    org = gh("api", f"orgs/{ORG}", "--jq", "{public_repos, followers, created_at}")
    counts = {}
    for name in ("org_participation.csv", "commit_weeks.csv", "repo_set.csv"):
        path = os.path.join(outdir, name)
        if not os.path.exists(path):
            continue
        frame = pd.read_csv(path)
        entry = {"rows": len(frame)}
        if "login" in frame.columns:
            entry["distinct_logins"] = int(frame["login"].nunique())
        if "created_at" in frame.columns and len(frame):
            entry["date_range"] = [str(frame["created_at"].min()),
                                   str(frame["created_at"].max())]
        counts[name] = entry
    snapshot = {"snapshot_date": date.today().isoformat(),
                "org": ORG, "main_repo": MAIN_REPO,
                "repo": repo, "org_meta": org, "files": counts}
    path = os.path.join(outdir, "snapshot.json")
    with open(path, "w") as f:
        json.dump(snapshot, f, indent=2)
    print(f"\n  snapshot {snapshot['snapshot_date']}:")
    for k, v in counts.items():
        extra = (f", {v['distinct_logins']} logins" if "distinct_logins" in v else "")
        span = (f", {v['date_range'][0]}..{v['date_range'][1]}" if "date_range" in v else "")
        print(f"    {k}: {v['rows']} rows{extra}{span}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--refresh", action="store_true",
                        help="re-fetch files that already exist")
    args = parser.parse_args()

    gh("api", "user", "--jq", "{login}")   # fails loudly if gh is not authenticated
    print(f"Collecting {ORG} community stats into data/raw/github_stats/")
    fetch_org_participation(args.refresh)
    repos = fetch_repo_set(args.refresh)["repo"].tolist()
    fetch_commit_weeks(repos, args.refresh)
    write_snapshot()


if __name__ == "__main__":
    main()
