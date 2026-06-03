# exomiser-rerank-benchmark

Does **OpenScientist's evidence-based reranking** improve on **Exomiser's** baseline
ranking for rare-disease gene prioritization?

This is a standalone evaluation harness. It is a **client of the OpenScientist REST API** —
it submits phenopackets as discovery jobs, then compares two rankings per case:

- **A — baseline:** Exomiser's `exomiser_ranking.json` (deterministic, tool-emitted).
- **B — reranked:** the agent's evidence-based reranking (PubMed / ClinVar / OMIM / db skills).

## Question & metric

> For a patient's HPO phenotypes, does reranking move the **true causal gene** higher than
> Exomiser's raw phenotype-driven ranking?

- **Mode:** Exomiser **phenotype-only** (no VCF — rank candidate genes by HPO match = the
  differential diagnosis).
- **Unit:** **gene** (true causal gene per case).
- **Primary metric:** **hits@k** for k ∈ {1, 3, 5, 10} — fraction of cases where the true gene
  is in the top-k — reported for baseline vs reranked.
- **Stats:** paired design → **McNemar's test** per k; plus **MRR** / mean-rank and per-case
  *helped / hurt / unchanged*.

## Data

Phenopackets with a known single causal gene, from **Monarch Phenopacket Store**
(<https://github.com/monarch-initiative/phenopacket-store>).

Two non-obvious requirements the harness enforces:

1. **Strip the diagnosis before submitting.** Phenopacket-Store files contain the answer
   (the `interpretations`/`diagnosis` block). Exomiser ignores it (it reads only HPO terms),
   but the LLM agent would *read it and cheat*. The harness feeds OS only
   `id + subject + phenotypicFeatures + metaData`, and keeps the ground-truth gene separately.
2. **Match genes by stable ID** (Entrez / Ensembl / HGNC), not symbol — symbols have aliases
   and would cause false misses.

## Harness flow

```
prepare.py  fetch phenopackets → extract (HPO, true gene id) → strip diagnosis → sanitized/ + ground_truth.csv
run.py      for each: POST /api/v1/jobs (sanitized phenopacket + prompt, max_iterations≈2) → poll → pull
            exomiser_ranking.json + the rerank output → results/<id>/
score.py    true-gene rank in baseline vs rerank → hits@k, McNemar, MRR → report + per-case CSV + plot
```

The reranking prompt asks the agent to run Exomiser (phenotype-only), then re-rank the top
candidate genes citing evidence, and emit a **structured** gene list (rank, geneSymbol, geneId,
exomiserRank, rationale, evidence) — the schema is defined in the prompt because OS treats
reranking as opt-in/unschematized.

## Watch-outs

- **Ceiling effect (run the pilot!):** clean phenopackets → Exomiser may already hit @1 → no
  headroom for reranking. The **10-case pilot** measures baseline hits@1 first; if it's ~90%+
  we focus on baseline-misses or add HPO noise before the full ~100.
- **LLM prior knowledge:** the model may "know" famous gene–disease links → inflated,
  non-generalizable gains. Mitigate by requiring citations + favoring rarer/recent diseases.
- **Nondeterminism:** reranking varies run-to-run; consider repeating a subset 3×.
- **Cost/where:** ~100 LLM jobs is hours of compute. Pilot runs locally against an
  exomiser-enabled OS (`:8080`, slow/emulated); the full run wants native compute (e.g. an
  OS instance on a server with Exomiser data installed).

## Setup

```bash
cp .env.example .env   # set OS_API_URL + OS_API_KEY (from the OS "API Keys" page)
uv sync
uv run python -m benchmark.prepare --n 10        # pilot: 10 cases
uv run python -m benchmark.run                   # submit + collect
uv run python -m benchmark.score                 # hits@k baseline vs rerank
```

`OS_API_KEY` and `.env` are gitignored — never commit them.

## Status

Scaffold. Module logic (`prepare` / `run` / `score`) is stubbed; see each module's docstring.
