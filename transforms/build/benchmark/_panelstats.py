"""Statistics for the common-panel validation — stratified AUC and its uncertainty.

The question is G1: do ECSPr edges light up with the conditions expected to move
them? Operationally: within a panel COLUMN (one product), does the condition that
is supposed to move it rank above the conditions that are not?

Why not BH over per-cell p-values
----------------------------------
It is arithmetically impossible, not merely unfashionable. The per-cell null here
has ~48 usable members after the degenerate-row gate, so the smallest attainable
empirical p is 1/49 = 0.0204. Bonferroni/BH over ~700 cells needs 0.05/700 =
7.1e-5. No cell can ever reject, whatever the biology. Reporting per-cell stars
off this null would be reporting the null's size.

So the unit of evidence is the COLUMN (a ranking), not the cell, and the estimand
is an AUC.

Why van Elteren rather than one pooled AUC
-------------------------------------------
Columns are ragged (per-facet liveness differs) and some share a positive (15 of
73 GOF columns, max 5). Pooling cells across columns would let a big column
dominate and would double-count shared positives. The stratified (van Elteren)
estimator handles both in one formula:

    AUC_S = Σ_j U_j / Σ_j |P_j|·|N_j|

— each stratum j contributes its own concordant-pair count and its own
normaliser, so a stratum's weight is its number of comparisons, and a stratum
with no positives or no negatives contributes nothing rather than a NaN.

Why permutation + cluster bootstrap rather than a closed form
--------------------------------------------------------------
The closed-form MW variance assumes independent observations. Cells are not
independent: one condition contributes a cell to EVERY column, so the condition —
not the cell — is the independent unit. The cluster bootstrap resamples
conditions; the permutation shuffles the condition→product assignment within a
facet, which is exactly the null "ECSPr does not know which product belongs to
which condition."

Env: any numpy. No graph, no solver — this module is deliberately testable alone
(`mamba run -n ml python _panelstats.py --selftest`).
"""
from __future__ import annotations

import numpy as np

# --- Mann-Whitney U with the tie term ---------------------------------------


def mannwhitney_u(pos: np.ndarray, neg: np.ndarray) -> float:
    """U = #{(p,n): p > n} + 0.5·#{(p,n): p == n}.

    The tie term is not optional here. After the degenerate-row gate the null
    still carries exact zeros (a condition can touch >=2 metabolites and still
    fail to move a given pair), and a route-creation column is all-zero by
    construction. Dropping ties would score those as wins.
    """
    pos = np.asarray(pos, dtype=np.float64)
    neg = np.asarray(neg, dtype=np.float64)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    gt = (pos[:, None] > neg[None, :]).sum()
    eq = (pos[:, None] == neg[None, :]).sum()
    return float(gt + 0.5 * eq)


def auc(pos: np.ndarray, neg: np.ndarray) -> float:
    """U / (|P|·|N|) — the probability a random positive outranks a random
    negative, ties counted as half."""
    pos, neg = np.asarray(pos), np.asarray(neg)
    if pos.size == 0 or neg.size == 0:
        return float("nan")
    return mannwhitney_u(pos, neg) / (pos.size * neg.size)


# --- stratified (van Elteren) AUC -------------------------------------------


def stratified_auc(strata: list) -> float:
    """AUC_S = Σ U_j / Σ |P_j|·|N_j| over strata [(pos_j, neg_j), ...].

    Strata with an empty side contribute nothing (no comparisons exist), rather
    than a NaN that would poison the sum.
    """
    num = 0.0
    den = 0.0
    for pos, neg in strata:
        pos, neg = np.asarray(pos), np.asarray(neg)
        if pos.size == 0 or neg.size == 0:
            continue
        num += mannwhitney_u(pos, neg)
        den += pos.size * neg.size
    return num / den if den > 0 else float("nan")


# --- percentile of a value against a null ------------------------------------


def null_percentile(x: float, null: np.ndarray) -> float:
    """(#{null <= x} + 0.5·#{null == x}) / |null| — the mid-rank percentile.

    Mid-ranks, not `<`, for the same reason `mannwhitney_u` counts ties: an
    all-zero null with x == 0 must score 0.5 ("indistinguishable"), not 1.0
    ("beat everything").
    """
    null = np.asarray(null, dtype=np.float64)
    if null.size == 0:
        return float("nan")
    lt = (null < x).sum()
    eq = (null == x).sum()
    return float(lt + 0.5 * eq) / null.size


def emp_p(x: float, null: np.ndarray) -> float:
    """(#{null >= x} + 1) / (|null| + 1) — the add-one empirical p.

    Floored at 1/(|null|+1) BY CONSTRUCTION. That floor is why these are
    reported as effect sizes and ranks, never thresholded — see the module
    docstring. Provided for diagnostics only.
    """
    null = np.asarray(null, dtype=np.float64)
    return float((null >= x).sum() + 1) / (null.size + 1)


# --- uncertainty --------------------------------------------------------------


def permutation_test(pos_by_stratum: dict, neg_pool: dict, n_perm: int = 10000,
                     seed: int = 0) -> tuple:
    """Permute which condition is the positive within each stratum's own pool.

    `pos_by_stratum`: {stratum_key: index of the true positive within pool}
    `neg_pool`:       {stratum_key: 1-D array of ALL conditions' values, positive
                       included at the index above}

    Returns (auc_obs, p_perm, null_aucs). The permutation reassigns the positive
    label within the stratum's own value pool — it does NOT resample values — so
    it holds the panel's value distribution fixed and tests only the
    condition→product assignment, which is the claim.
    """
    rng = np.random.default_rng(seed)
    keys = sorted(pos_by_stratum)

    def _auc_for(assign: dict) -> float:
        strata = []
        for k in keys:
            vals = np.asarray(neg_pool[k], dtype=np.float64)
            i = assign[k]
            mask = np.ones(vals.size, dtype=bool)
            mask[i] = False
            strata.append((vals[i:i + 1], vals[mask]))
        return stratified_auc(strata)

    auc_obs = _auc_for(pos_by_stratum)
    null = np.empty(n_perm)
    for b in range(n_perm):
        assign = {k: int(rng.integers(0, len(neg_pool[k]))) for k in keys}
        null[b] = _auc_for(assign)
    p = float((null >= auc_obs).sum() + 1) / (n_perm + 1)
    return auc_obs, p, null


def cluster_bootstrap(strata_by_condition: dict, n_boot: int = 2000,
                      seed: int = 0) -> tuple:
    """Resample CONDITIONS (with replacement) and recompute AUC_S.

    `strata_by_condition`: {condition_id: [(pos_j, neg_j), ...]} — the strata a
    condition contributes. The condition is the independent unit: it appears in
    every column, so resampling cells would understate the interval.

    Returns (lo, hi) percentile-95 bounds.
    """
    rng = np.random.default_rng(seed)
    conds = sorted(strata_by_condition)
    out = np.empty(n_boot)
    for b in range(n_boot):
        pick = rng.integers(0, len(conds), size=len(conds))
        strata = []
        for i in pick:
            strata.extend(strata_by_condition[conds[i]])
        out[b] = stratified_auc(strata)
    out = out[np.isfinite(out)]
    if out.size == 0:
        return float("nan"), float("nan")
    return float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5))


# --- self-test ---------------------------------------------------------------


def _selftest() -> int:
    ok = True

    def check(name, got, want, tol=1e-12):
        nonlocal ok
        good = (np.isnan(got) and np.isnan(want)) or abs(got - want) <= tol
        print(f"  {'ok  ' if good else 'FAIL'} {name}: got {got!r} want {want!r}")
        if not good:
            ok = False

    print("mannwhitney_u / auc — known answers")
    # perfect separation
    check("auc perfect", auc([3, 4, 5], [0, 1, 2]), 1.0)
    check("auc reversed", auc([0, 1, 2], [3, 4, 5]), 0.0)
    # all ties -> 0.5, the tie term's whole point
    check("auc all ties", auc([0, 0, 0], [0, 0, 0]), 0.5)
    # a route-creation column: positive 0, null all 0 -> must be 0.5, not 1.0
    check("auc zero-vs-zero", auc([0.0], [0.0] * 48), 0.5)
    # half ties
    check("auc half tie", auc([1.0], [0.0, 1.0]), 0.75)
    # empty side
    check("auc empty pos", auc([], [1, 2]), float("nan"))

    print("stratified_auc")
    # two identical strata == one of them
    s = [([3.0], [1.0, 2.0]), ([3.0], [1.0, 2.0])]
    check("two identical strata", stratified_auc(s), 1.0)
    # ragged: stratum A perfect (1 vs 2), stratum B reversed (1 vs 4).
    # U_A=2, U_B=0 -> 2/(2+4) = 0.3333...
    s = [([3.0], [1.0, 2.0]), ([0.0], [1.0, 2.0, 3.0, 4.0])]
    check("ragged weighting", stratified_auc(s), 2.0 / 6.0)
    # an empty stratum contributes nothing rather than NaN
    s = [([3.0], [1.0, 2.0]), ([], [1.0])]
    check("empty stratum ignored", stratified_auc(s), 1.0)
    check("all strata empty", stratified_auc([([], [])]), float("nan"))

    print("null_percentile / emp_p")
    check("pct above all", null_percentile(10, [1, 2, 3]), 1.0)
    check("pct below all", null_percentile(0, [1, 2, 3]), 0.0)
    check("pct midrank tie", null_percentile(2, [1, 2, 3]), 0.5)
    check("pct zero vs zeros", null_percentile(0.0, [0.0] * 10), 0.5)
    # the floor that makes BH impossible — stated as a test so it cannot be forgotten
    check("emp_p floor at 1/(N+1)", emp_p(1e9, np.zeros(48)), 1.0 / 49.0)

    print("cluster_bootstrap — degenerate input is honest")
    lo, hi = cluster_bootstrap({"c1": [([3.0], [1.0, 2.0])]}, n_boot=100)
    check("bootstrap of a constant", lo, 1.0)
    check("bootstrap of a constant (hi)", hi, 1.0)

    print("permutation_test — CALIBRATION on noise (not a single panel)")
    # Asserting that ONE noise panel is non-significant is not a test of the
    # method; it is a test of the seed. With 12 strata, AUC_S of a noise panel
    # has sd ~ 0.29/sqrt(12) ~ 0.083, so ~1 noise panel in 20 lands beyond 2 sd
    # and SHOULD return a small p — a test that forbids that is demanding a
    # miscalibrated test. (Measured: seed 1 gives AUC=0.730, p=0.005, and it is
    # a correct answer for that draw.) The property that must hold is that p is
    # Uniform(0,1) under the null.
    rng = np.random.default_rng(0)
    n_panels = 120
    ps = np.empty(n_panels)
    for t in range(n_panels):
        pool = {f"s{j}": rng.normal(size=30) for j in range(12)}
        _, p, _ = permutation_test({k: 0 for k in pool}, pool, n_perm=200, seed=t)
        ps[t] = p
    for alpha in (0.05, 0.25, 0.50):
        frac = float((ps <= alpha).mean())
        se = np.sqrt(alpha * (1 - alpha) / n_panels)
        good = abs(frac - alpha) < 4 * se           # 4 se: a gate, not a coin flip
        print(f"  {'ok  ' if good else 'FAIL'} P(p<={alpha:.2f}) = {frac:.3f} "
              f"(expect {alpha:.2f} +/- {4*se:.3f})")
        if not good:
            ok = False

    print("permutation_test — a planted signal must be significant")
    pool2 = {}
    for j in range(12):
        v = rng.normal(size=30)
        v[0] = 10.0                                   # plant the positive on top
        pool2[f"s{j}"] = v
    a2, p2, _ = permutation_test({k: 0 for k in pool2}, pool2, n_perm=2000, seed=3)
    print(f"  planted panel: AUC={a2:.3f} p={p2:.3f}")
    if a2 < 0.99 or p2 > 0.001:
        print("  FAIL: planted signal not recovered")
        ok = False
    else:
        print("  ok   planted signal recovered")

    print()
    print("SELFTEST PASS" if ok else "SELFTEST FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    import sys
    if "--selftest" in sys.argv:
        sys.exit(_selftest())
    print(__doc__)
