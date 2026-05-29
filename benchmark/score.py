"""Score baseline (Exomiser) vs reranked (OS) against the true causal gene.

For each case: the true gene's 1-based rank in each ranking (symbol match; HGNC kept in
ground_truth for a future stable-ID mapping if alias misses appear). Then hits@{1,3,5,10}
for both conditions, McNemar's paired exact test per k, and MRR.

Outputs a table to stdout and results/scores.csv.
"""

from __future__ import annotations

import csv
import json

from statsmodels.stats.contingency_tables import mcnemar

from . import config

K_VALUES = (1, 3, 5, 10)


def _norm(s: object) -> str:
    return str(s or "").strip().upper()


def _load_ranking(case_id: str, name: str) -> list[dict] | None:
    p = config.RESULTS / case_id / name
    if not p.exists():
        return None
    data = json.loads(p.read_text())
    if name == "exomiser_ranking.json":
        return data.get("geneRanking", [])  # tool schema: {geneRanking: [...], ...}
    if isinstance(data, list):  # reranked.json may be a bare array
        return data
    return data.get("ranking", [])  # ...or {ranking: [...]}


def _rank_of(true_symbol: str, ranking: list[dict] | None) -> int | None:
    if ranking is None:
        return None
    target = _norm(true_symbol)
    for i, entry in enumerate(ranking, 1):
        if _norm(entry.get("geneSymbol")) == target:
            return i
    return None


def main() -> None:
    if not config.GROUND_TRUTH.exists():
        raise SystemExit("ground_truth.csv not found")
    cases = list(csv.DictReader(config.GROUND_TRUTH.open()))

    per_case = []
    for c in cases:
        cid, sym = c["case_id"], c["gene_symbol"]
        br = _rank_of(sym, _load_ranking(cid, "exomiser_ranking.json"))
        rr = _rank_of(sym, _load_ranking(cid, "reranked.json"))
        per_case.append(
            {
                "case_id": cid,
                "gene": sym,
                "baseline_rank": br,
                "reranked_rank": rr,
                "outcome": _outcome(br, rr),
            }
        )

    n = len(per_case)
    print(f"\ncases scored: {n}\n")
    print(f"{'k':>3} | {'baseline':>12} | {'reranked':>12} | {'McNemar p':>10}")
    print("-" * 50)
    for k in K_VALUES:
        b = [
            pc["baseline_rank"] is not None and pc["baseline_rank"] <= k
            for pc in per_case
        ]
        r = [
            pc["reranked_rank"] is not None and pc["reranked_rank"] <= k
            for pc in per_case
        ]
        both = sum(1 for x, y in zip(b, r) if x and y)
        b_only = sum(1 for x, y in zip(b, r) if x and not y)
        r_only = sum(1 for x, y in zip(b, r) if not x and y)
        neither = sum(1 for x, y in zip(b, r) if not x and not y)
        try:
            p = mcnemar([[both, b_only], [r_only, neither]], exact=True).pvalue
        except (ValueError, ZeroDivisionError):
            p = float("nan")
        print(
            f"{k:>3} | {sum(b):>4}/{n} ({sum(b) / n:>4.0%}) | "
            f"{sum(r):>4}/{n} ({sum(r) / n:>4.0%}) | {p:>10.3f}"
        )

    def _mrr(key: str) -> float:
        return sum((1.0 / pc[key]) if pc[key] else 0.0 for pc in per_case) / n

    print(
        f"\nMRR: baseline={_mrr('baseline_rank'):.3f}  reranked={_mrr('reranked_rank'):.3f}"
    )
    moved = sum(pc["outcome"] == "helped" for pc in per_case)
    hurt = sum(pc["outcome"] == "hurt" for pc in per_case)
    print(
        f"reranking moved true gene UP in {moved}, DOWN in {hurt}, else unchanged/absent\n"
    )

    out = config.RESULTS / "scores.csv"
    config.RESULTS.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(per_case[0].keys()))
        writer.writeheader()
        writer.writerows(per_case)
    print(f"per-case -> {out}")


def _outcome(br: int | None, rr: int | None) -> str:
    if br is None and rr is None:
        return "absent"
    if rr is not None and (br is None or rr < br):
        return "helped"
    if br is not None and (rr is None or rr > br):
        return "hurt"
    return "unchanged"


if __name__ == "__main__":
    main()
