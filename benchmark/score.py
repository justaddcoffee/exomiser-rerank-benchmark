"""Score baseline vs reranked against ground truth and run the stats.

For each case: find the true causal gene's rank in (a) Exomiser's geneRanking and (b) the
reranked list, matching on stable gene ID (Entrez/Ensembl/HGNC), falling back to symbol.

Outputs:
  - hits@k for k in {1,3,5,10}, baseline vs reranked, with McNemar's paired test per k.
  - MRR / mean-rank for each condition.
  - per-case CSV: id, true_gene, baseline_rank, reranked_rank, delta, helped/hurt/unchanged.

TODO(impl):
  - rank_of(true_gene_ids, ranking) -> int | None  (stable-ID match, symbol fallback)
  - hits_at_k, mrr
  - mcnemar paired test (statsmodels.stats.contingency_tables.mcnemar) per k
  - emit results table + per-case CSV (+ optional bar plot)
"""

from __future__ import annotations

K_VALUES = (1, 3, 5, 10)


def main() -> None:
    raise NotImplementedError("score: hits@k baseline vs reranked + McNemar over results/")


if __name__ == "__main__":
    main()
