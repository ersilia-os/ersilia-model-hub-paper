"""Additively weighted Voronoi (power) diagrams and area-controlled Voronoi treemaps.

Generic geometry, no plotting and no new dependencies — the plot classes live in
:mod:`plots_metadata`. Implements the Balzer & Deussen Voronoi-treemap scheme: iterate a
power diagram, moving each site to its cell centroid (Lloyd relaxation) and adjusting its
weight until each cell's area matches a target share of the enclosing polygon.

A *power* diagram is used rather than a multiplicatively weighted one because its bisectors
are straight lines: the cell of site ``i`` is
``{x : |x - p_i|^2 - w_i <= |x - p_j|^2 - w_j for all j}``, and expanding that inequality
cancels the quadratic term, leaving a half-plane. Every cell is therefore a convex polygon
obtainable by clipping the boundary with one half-plane per other site — exact, and directly
usable as a vector path.

Weights are fitted by gradient ascent on the semi-discrete optimal-transport dual, which is
concave with gradient ``goal_i - area_i``. That replaces the paper's multiplicative update on
non-negative "radii" with an overlap clamp: the clamp shrinks every weight toward parity, so
targets more than a few-fold apart never separate, whereas ascent on the dual lets weights go
negative and reaches <1% area error even at a 1:100,000 ratio. See ``voronoi_treemap`` for the
two details that make it work — preconditioning the gradient by ``1/goal``, and not treating a
failed line search as convergence.
"""

import numpy as np


def polygon_area(poly):
    """Shoelace area of a closed-by-implication polygon (>= 0); 0 for degenerate input."""
    if len(poly) < 3:
        return 0.0
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))


def polygon_centroid(poly):
    """Area centroid of a polygon; falls back to the vertex mean when degenerate."""
    if len(poly) < 3:
        return poly.mean(axis=0)
    x, y = poly[:, 0], poly[:, 1]
    cross = x * np.roll(y, -1) - np.roll(x, -1) * y
    a = cross.sum()
    if abs(a) < 1e-15:
        return poly.mean(axis=0)
    cx = ((x + np.roll(x, -1)) * cross).sum() / (3.0 * a)
    cy = ((y + np.roll(y, -1)) * cross).sum() / (3.0 * a)
    return np.array([cx, cy])


def clip_halfplane(poly, a, b, c):
    """Sutherland-Hodgman clip of ``poly`` to the half-plane ``a*x + b*y <= c``."""
    if len(poly) == 0:
        return poly
    d = a * poly[:, 0] + b * poly[:, 1] - c        # <= 0 is inside
    out = []
    n = len(poly)
    for i in range(n):
        j = (i + 1) % n
        di, dj = d[i], d[j]
        if di <= 0:
            out.append(poly[i])
        if (di < 0 < dj) or (dj < 0 < di):          # edge crosses the boundary
            t = di / (di - dj)
            out.append(poly[i] + t * (poly[j] - poly[i]))
    return np.array(out) if out else np.empty((0, 2))


def power_cells(sites, weights, boundary):
    """Power-diagram cells of ``sites`` (weights ``w``) clipped to the ``boundary`` polygon.

    Returns one array of vertices per site, in site order; an array with fewer than 3 rows
    means that site was swallowed by its neighbours and has no cell.
    """
    n = len(sites)
    sq = (sites ** 2).sum(axis=1)
    cells = []
    for i in range(n):
        poly = boundary
        for j in range(n):
            if i == j or len(poly) < 3:
                continue
            # |x-pi|^2 - wi <= |x-pj|^2 - wj  ->  2(pj-pi).x <= |pj|^2 - wj - |pi|^2 + wi
            a, b = 2.0 * (sites[j, 0] - sites[i, 0]), 2.0 * (sites[j, 1] - sites[i, 1])
            c = sq[j] - weights[j] - sq[i] + weights[i]
            poly = clip_halfplane(poly, a, b, c)
        cells.append(poly)
    return cells


def quadratic_moment(poly, p):
    """``integral over poly of |x - p|^2 dx``, exactly.

    The polygon is fanned into triangles from its first vertex (valid because power cells
    clipped from a convex boundary are convex), and on a triangle the three-edge-midpoint rule
    is exact for a quadratic integrand.
    """
    if len(poly) < 3:
        return 0.0
    total = 0.0
    v0 = poly[0]
    for k in range(1, len(poly) - 1):
        a, b, c = v0, poly[k], poly[k + 1]
        area = 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))
        if area == 0.0:
            continue
        mids = 0.5 * np.array([a + b, b + c, c + a])
        total += area / 3.0 * ((mids - p) ** 2).sum()
    return total


def _seed_sites(boundary, n, rng):
    """``n`` random points inside ``boundary`` (rejection sampling on its bounding box)."""
    lo, hi = boundary.min(axis=0), boundary.max(axis=0)
    pts = []
    while len(pts) < n:
        for p in rng.uniform(lo, hi, size=(4 * n + 16, 2)):
            if len(pts) == n:
                break
            if _inside(boundary, p):
                pts.append(p)
    return np.array(pts[:n], dtype=float)


def _inside(poly, p):
    """Ray-casting point-in-polygon test."""
    x, y = p
    xs, ys = poly[:, 0], poly[:, 1]
    xj, yj = np.roll(xs, -1), np.roll(ys, -1)
    crosses = ((ys > y) != (yj > y)) & (x < xs + (y - ys) / (yj - ys + 1e-300) * (xj - xs))
    return bool(crosses.sum() % 2)


def voronoi_treemap(targets, boundary, *, seed, iterations=1200, tol=0.01, settle=0.6,
                    restarts=4):
    """Area-controlled Voronoi treemap, retried from several seeds until it fits.

    Thin wrapper over :func:`_fit_treemap`: convergence depends on where the sites happen to
    start (a 1:11,000 spread over 12 cells fails from some seeds and fits from others), so this
    tries up to ``restarts`` seeds and keeps the best fit, stopping as soon as one meets ``tol``.
    Returns ``(cells, error)``; see :func:`_fit_treemap` for the arguments and for why the error
    must be checked rather than assumed small.
    """
    best_cells, best_err = None, np.inf
    for k in range(max(1, restarts)):
        cells, err = _fit_treemap(targets, boundary, seed=seed + 1000 * k,
                                  iterations=iterations, tol=tol, settle=settle)
        if err < best_err:
            best_cells, best_err = cells, err
        if best_err < tol:
            break
    return best_cells, best_err


def _fit_treemap(targets, boundary, *, seed, iterations=1200, tol=0.01, settle=0.6):
    """Tessellate ``boundary`` into one convex cell per target, area proportional to it.

    For fixed sites, finding weights that give prescribed cell areas is a concave maximisation
    whose gradient in ``w_i`` is exactly ``goal_i - area_i`` (the semi-discrete optimal-transport
    dual), so the weights are fitted by gradient ascent with a backtracking step. Unlike a
    multiplicative update on non-negative "radii", this needs no overlap clamp — a swallowed
    cell has area 0, hence a positive gradient, and grows back on its own — and it lets weights
    go negative, which is what allows cells orders of magnitude apart in size.

    Site positions are relaxed toward their cell centroids (Lloyd) for the first ``settle``
    fraction of the iterations, which is what makes the cells compact; positions then freeze so
    the remaining iterations converge the areas against a fixed problem.

    Parameters
    ----------
    targets    : positive weights; each cell's area converges to ``target / sum(targets)``
                 of the boundary's area.
    boundary   : ``(k, 2)`` array of polygon vertices. Cells are convex only if it is.
    seed       : RNG seed for the initial site positions — the one stochastic input, so the
                 layout is reproducible.
    iterations : maximum rounds.
    tol        : stop once the largest relative area error falls below this.
    settle     : fraction of iterations during which sites are still allowed to move.

    Returns ``(cells, error)`` — the polygon per target in input order, and the largest
    relative area error ``max(|area - goal| / goal)`` achieved. **Check the error**: targets
    spanning many orders of magnitude cannot all be honoured, since a cell far below the site
    spacing cannot be resolved, and the caller should report that rather than imply the areas
    are exact.
    """
    targets = np.asarray(targets, dtype=float)
    n = len(targets)
    total = polygon_area(boundary)
    goal = targets / targets.sum() * total

    if n == 1:
        return [np.asarray(boundary, dtype=float)], 0.0

    rng = np.random.default_rng(seed)
    sites = _seed_sites(boundary, n, rng)
    # Weights behave like squared radii, so seeding them with the goal areas starts the fit near
    # the right spread. Starting from zeros is an unweighted diagram — every cell roughly equal —
    # which is a poor place to begin when the targets span orders of magnitude.
    w = goal / np.pi
    w -= w.mean()
    step = 0.5 * total                # weights have units of area, so the step scale does too
    stalled = 0
    step_cycle = [0.5, 0.05, 2.0, 0.005]

    def evaluate(weights):
        """Cells, their areas, the concave dual objective, and the max relative area error.

        The objective maximised is ``sum_i [ int_{C_i} (|x-p_i|^2 - w_i) dx + w_i * goal_i ]``,
        whose gradient in ``w_i`` is ``goal_i - area_i`` by the envelope theorem. It is the
        line-search objective precisely because the max relative error is not usable for one:
        the moment a cell is swallowed that error pins to 1.0 and no step can improve it, so a
        search on it halves the step to nothing and the fit stalls with the cell still empty.
        """
        cells = power_cells(sites, weights, boundary)
        areas = np.empty(n)
        obj = 0.0
        for i, cell in enumerate(cells):
            areas[i] = polygon_area(cell)
            obj += quadratic_moment(cell, sites[i]) + weights[i] * (goal[i] - areas[i])
        return cells, areas, obj, float(np.max(np.abs(areas - goal) / goal))

    cells, areas, obj, err = evaluate(w)
    best_cells, best_err = cells, err

    for it in range(iterations):
        if err < tol:
            break

        # Preconditioned ascent: the raw gradient is goal_i - area_i, so a cell whose target is
        # 1000x smaller moves 1000x slower and never catches up. Dividing through by goal_i
        # turns it into a relative area error, equalising the rates. Scaling by any positive
        # diagonal preserves the ascent direction on a concave objective, so the line search
        # below is still valid.
        grad = (goal - areas) / goal
        grad = grad / max(np.abs(grad).max(), 1e-300)   # unit max, so `step` means one thing
        advanced = False
        for _ in range(24):
            trial = w + step * grad
            trial -= trial.mean()           # only weight differences matter; keep them centred
            t_cells, t_areas, t_obj, t_err = evaluate(trial)
            if t_obj > obj:
                w, cells, areas, obj, err = trial, t_cells, t_areas, t_obj, t_err
                step = min(step * 1.5, 4.0 * total)
                advanced = True
                break
            step *= 0.4
        if not advanced:
            # A failed line search means the step is out of scale, not that we are done. Cycle
            # through a few scales before giving up, since the usable range varies a lot with
            # how unequal the targets are.
            step = step_cycle[stalled % len(step_cycle)] * total
            stalled += 1
            if stalled >= 3 * len(step_cycle):
                break
        else:
            stalled = 0

        if err < best_err:
            best_err, best_cells = err, cells

        if it < settle * iterations:
            moved = sites.copy()
            spare = int(np.argmax(areas - goal))     # the most over-sized cell
            empty = [i for i, a in enumerate(areas) if a <= 0]
            for i, cell in enumerate(cells):
                if areas[i] > 0:
                    moved[i] = polygon_centroid(cell)
            if empty and len(cells[spare]) >= 3:
                # Swallowed cells are re-dropped into the cell with the most area to give away,
                # at DISTINCT random points. Sending them all to its centroid instead makes them
                # effectively coincident, and coincident sites just swallow one another — which
                # is what left the 12-cell, 1:11,000 E. coli case stuck with 9 empty cells.
                fresh = _seed_sites(cells[spare], len(empty), rng)
                for i, p in zip(empty, fresh):
                    moved[i] = p
            sites = moved
            cells, areas, obj, err = evaluate(w)
            if err < best_err:
                best_err, best_cells = err, cells

    return best_cells, float(best_err)
