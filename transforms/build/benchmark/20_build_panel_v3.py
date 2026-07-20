"""Emit v3's panel -- Y/panel_edges.tsv and Y/panel_anchors.tsv -- and gate every anchor.

The panel is the benchmark's coordinate system: every condition is measured
against ONE common set of `anchor x product` edges, so a column is a ranking and
the claim becomes "the right condition tops it".

WHY THIS SCRIPT HARD-FAILS RATHER THAN WARNS
--------------------------------------------
A dead anchor is SILENT. Conductance between two metabolites is 0.0 when an
endpoint is off the largest connected component -- the same value it returns for
"connected, no conductance" -- so an anchor with a wrong MetaNetX id loses every
condition beneath it and reports a column of honest-looking zeros. Three ids in
this project were already dead this way (G6P, G3P, sulfide). Liveness is a GATE,
not a diagnostic, and every id is re-measured on every facet.

WHAT LIVENESS IS *NOT*
----------------------
Anchor liveness gates the PANEL, which is X-side coordinate machinery. It must
never gate an EXPECTATION. Y declares biology; whether a network can represent
cardiolipin is the implementation's problem, reported by the scorer as a coverage
miss. The moment a Y row is dropped because the incumbent cannot score it, Y has
been tuned to the incumbent. So this emits `live_<facet>` columns and lets the
scorer report coverage -- it does not delete edges.

WHY THE PANEL EXTENDS v1's RATHER THAN RE-DERIVING IT
-----------------------------------------------------
v1's `products()` unions four sources, two of which (`_gof_specs.GOF_SPECS` and
`_panelspec`) live in the metabolic-modelling scope. Reaching into them would
break the self-containment that is the entire point of the v3 tree.

So the product set is v1's FROZEN product set -- staged here as
`contract/panel_edges.tsv`, which already encodes all four of those sources --
UNION v3's own new targets. This is a legitimate extension, not an inheritance
of conclusions: the panel is coordinate machinery, and `contract/` is declared
as exactly that, the contract SHAPE. No EXPECTATION is read from v1; those are
re-derived from ground truth in 21_build_Y_v3.py.

The anchors are treated the same way -- v1's 49 ids are read as ids, then
RE-GATED live against v3's own declared base graphs. The gate is the load-bearing
part and it re-runs in full.

WHERE THE BASE GRAPHS COME FROM
-------------------------------
`canon.BENCH_V3_BASE_GRAPHS`, a declared library item. They are NOT rebuilt from
X, and cannot be: X withholds the atom mapping on purpose, and the atom-transit
count is both the edge weight and -- via the w<=0 edge drop -- the connectivity.
Base graphs are network construction, never ECSPr output, which is what lets this
builder read them without breaking answer-key independence.

Usage:  mamba run -n ml python 20_build_panel_v3.py
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import networkx as nx
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _v3bench import (  # noqa: E402
    CONTRACT_DIR, CURRENCY_NAMES, ELEMENTS, FACETS, OBS_DIR, X_DIR, Y_DIR,
    canon, write_tsv,
)

BASE_DIR = Path(canon.BENCH_V3_BASE_GRAPHS)
TARGETS = Path(canon.BENCH_V3_DECISIONS) / "target_resolution.tsv"


def load_base(facet: str, el: str) -> nx.Graph:
    with (BASE_DIR / facet / f"base_{el}.pkl").open("rb") as fh:
        return pickle.load(fh)


def lcc_of(facet: str, el: str) -> tuple[set, nx.Graph]:
    g = load_base(facet, el)
    lcc = max(nx.connected_components(g), key=len)
    return lcc, g.subgraph(lcc)


def element_map() -> dict[str, set[str]]:
    """mnxm -> the elements of {C,N,S,P} its formula actually contains.

    A product is only a panel product on an element it CARRIES. Registering
    sulfate as a carbon product would put an edge in the panel that can never
    respond, which reads downstream as a failure to respond rather than as a
    cell that was never meaningful.

    Parsed by scanning for the element symbol not followed by a lowercase
    letter, so the S in `Se` and the C in `Cl` and `Co` do not count.
    """
    import re
    met = pd.read_csv(X_DIR / "metabolites.tsv", sep="\t", dtype=str).fillna("")
    pats = {el: re.compile(el + r"(?![a-z])") for el in ELEMENTS}
    out = {}
    for mnxm, formula in zip(met["mnxm"], met["formula"]):
        if not formula:
            continue
        out[mnxm] = {el for el, p in pats.items() if p.search(formula)}
    return out


def currency_mnxm() -> set[str]:
    met = pd.read_csv(X_DIR / "metabolites.tsv", sep="\t", dtype=str).fillna("")
    return {m for m, n in zip(met["mnxm"], met["name"])
            if n.strip().lower() in CURRENCY_NAMES}


def v3_products(elmap: dict, in_x: set) -> pd.DataFrame:
    """v3's own targets: the mechanical near-products, plus the resolved names.

    TWO SOURCES, and the order matters because it is the order of authority:

      1. MECHANICAL -- the direct non-currency products of every reaction a GOF
         condition INSERTS. This is v1's own rule and it is what carries the arm:
         380 of 382 GOF rows have a perturbed reaction present in X. It needs no
         name matching at all.
      2. RESOLVED   -- the free-text `target` / `source_metabolite_raw` strings
         resolved to ids in target_resolution.tsv. This is ENRICHMENT. Its
         coverage is intrinsically partial (heterologous products outside X,
         compound classes with no single id) and that is reported, not hidden.
    """
    curr = currency_mnxm()
    rp = pd.read_csv(X_DIR / "reaction_participants.tsv", sep="\t", dtype=str).fillna("")
    prod_of = (rp[rp.role == "product"]
               .groupby("mnxr")["mnxm"].apply(list).to_dict())

    rows = []

    gof = pd.read_csv(OBS_DIR / "gof_observations.tsv", sep="\t", dtype=str).fillna("")
    for mnxr_cell in gof["add_mnxr"]:
        for mnxr in (m.strip() for m in mnxr_cell.split(";") if m.strip()):
            for mnxm in prod_of.get(mnxr, []):
                if mnxm in curr or mnxm not in in_x:
                    continue
                for el in elmap.get(mnxm, ()):
                    rows.append(dict(element=el, mnxm=mnxm, source="v3_gof_near",
                                     category="inserted_reaction_product"))

    tr = pd.read_csv(TARGETS, sep="\t", dtype=str).fillna("")
    for _, r in tr[tr.mnxm != ""].iterrows():
        if r.mnxm not in in_x:
            continue
        for el in elmap.get(r.mnxm, ()):
            rows.append(dict(element=el, mnxm=r.mnxm,
                             source=f"v3_{r.arm}_target", category=r.method))

    return pd.DataFrame(rows)


def main() -> int:
    Y_DIR.mkdir(parents=True, exist_ok=True)

    # ---- gate: every anchor live in every facet's LCC ------------------------
    anc = pd.read_csv(CONTRACT_DIR / "panel_anchors.tsv", sep="\t")
    keep = ["element", "anchor", "mnxm", "label", "tier", "carrier_contaminated"]
    anc = anc[keep]

    print("ANCHOR LIVENESS GATE\n" + "=" * 78)
    dead, gate_rows = [], []
    for el in ELEMENTS:
        subs = {f: lcc_of(f, el) for f in FACETS}
        n_met = ", ".join(
            f"{f.replace('net', '')} {sum(1 for n in subs[f][0] if n[0] == 'met')}"
            for f in FACETS)
        print(f"\n  element {el}  (LCC mets: {n_met})")
        for _, a in anc[anc.element == el].iterrows():
            row = a.to_dict()
            state = []
            for f in FACETS:
                lcc, sub = subs[f]
                node = ("met", a.mnxm)
                if node in lcc:
                    deg = sub.degree(node)
                    w = sum(d.get(f"w_{el}", 0) for _, _, d in sub.edges(node, data=True))
                    row[f"deg_{f}"], row[f"w_{f}"], row[f"live_{f}"] = deg, f"{w:g}", True
                    state.append(f"d{deg}")
                else:
                    row[f"deg_{f}"], row[f"w_{f}"], row[f"live_{f}"] = 0, "0", False
                    dead.append(f"{el}/{a.anchor} ({a.mnxm}) dead on {f}")
                    state.append("DEAD")
            gate_rows.append(row)
            flag = " [carrier]" if a.carrier_contaminated else ""
            print(f"    {a.tier:<9}{a.anchor:<13}{a.mnxm:<14}"
                  f"{' '.join(f'{s:>7}' for s in state)}{flag}")

    if dead:
        print("\n" + "!" * 78)
        for d in dead:
            print(f"  DEAD ANCHOR: {d}")
        print("A dead anchor reports a column of honest-looking zeros. Fix the id.")
        print("!" * 78)
        return 1
    print(f"\n  GATE PASSED -- {len(gate_rows)} anchors live on all {len(FACETS)} facets")
    write_tsv(pd.DataFrame(gate_rows), Y_DIR / "panel_anchors.tsv")

    # ---- products: v1's frozen set, extended by v3's own --------------------
    met = pd.read_csv(X_DIR / "metabolites.tsv", sep="\t", dtype=str).fillna("")
    in_x = set(met["mnxm"])
    names = dict(zip(met["mnxm"], met["name"]))
    elmap = element_map()

    v1e = pd.read_csv(CONTRACT_DIR / "panel_edges.tsv", sep="\t", dtype=str)
    carried = (v1e[["element", "product_mnxm"]].drop_duplicates()
               .rename(columns={"product_mnxm": "mnxm"}))
    carried["source"] = "v1_panel"
    carried["category"] = "carried"

    new = v3_products(elmap, in_x)
    allp = pd.concat([carried, new], ignore_index=True)
    prods = (allp.groupby(["element", "mnxm"])
                 .agg(sources=("source", lambda s: ",".join(sorted(set(s)))),
                      categories=("category", lambda s: ",".join(sorted(set(s)))))
                 .reset_index())
    prods["name"] = prods["mnxm"].map(names).fillna("")

    n_new = len(prods) - len(carried.drop_duplicates(["element", "mnxm"]))
    print(f"\nPRODUCTS\n{'=' * 78}")
    print(f"  carried from the v1 panel : {len(carried.drop_duplicates(['element','mnxm'])):>5}")
    print(f"  added by v3               : {n_new:>5}")
    print(f"  panel product set         : {len(prods):>5}")

    # ---- the panel ----------------------------------------------------------
    edges = []
    for el in ELEMENTS:
        pr = prods[prods.element == el]
        for _, a in anc[anc.element == el].iterrows():
            for _, p in pr.iterrows():
                if p.mnxm == a.mnxm:
                    # the scorer excludes anchor==product; such a cell is a
                    # structural blank, not a measurement.
                    continue
                edges.append(dict(
                    edge_id=f"{el}__{a.anchor}__{p.mnxm}",
                    element=el, anchor=a.anchor, anchor_mnxm=a.mnxm,
                    anchor_tier=a.tier, anchor_carrier=a.carrier_contaminated,
                    product_mnxm=p.mnxm, product_name=p["name"],
                    product_sources=p.sources, product_categories=p.categories,
                ))
    panel = pd.DataFrame(edges)

    # product liveness per facet -- a coverage OUTPUT, never a filter
    for el in ELEMENTS:
        for f in FACETS:
            lcc, _ = lcc_of(f, el)
            live = {n[1] for n in lcc if n[0] == "met"}
            m = panel.element == el
            panel.loc[m, f"product_live_{f}"] = panel.loc[m, "product_mnxm"].isin(live)

    print(f"\nPANEL\n{'=' * 78}")
    for el in ELEMENTS:
        sub = panel[panel.element == el]
        cov = [f"{f.replace('net', '')} {100 * sub[f'product_live_{f}'].mean():.0f}%"
               for f in FACETS]
        print(f"  {el}: {sub.anchor.nunique():>2} anchors x {sub.product_mnxm.nunique():>3}"
              f" products = {len(sub):>5} edges   product live: {' / '.join(cov)}")

    write_tsv(panel, Y_DIR / "panel_edges.tsv")

    (Y_DIR / "_panel_provenance.json").write_text(json.dumps(dict(
        version="v3",
        facets=FACETS,
        n_anchors=int(len(gate_rows)),
        n_edges=int(len(panel)),
        per_element={el: int((panel.element == el).sum()) for el in ELEMENTS},
        product_sources=[
            "v1 frozen panel product set (set4 axes + LOF auxotrophy + GOF near-products)",
            "v3 mechanical near-products: non-currency products of every inserted reaction",
            "v3 resolved free-text targets (target_resolution.tsv)",
        ],
        base_graphs=str(BASE_DIR),
        gate=("every anchor re-measured live in every facet's LCC against v3's own "
              "declared base graphs; exits 1 on a dead anchor"),
        liveness_policy=("product liveness is a coverage OUTPUT, never a filter -- "
                         "Y declares biology, the scorer reports representability"),
        independence=("base graphs are network construction, not ECSPr output; no "
                      "expectation is read from v1, only the product coordinate set"),
    ), indent=2) + "\n")

    print(f"\npanel frozen -- {len(panel):,} edges over {len(prods):,} products")
    return 0


if __name__ == "__main__":
    sys.exit(main())
