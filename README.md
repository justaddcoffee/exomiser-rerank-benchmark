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

Two **corpora**, selected throughout the harness with `--corpus` (default `store`):

- **`store`** (pilot 1): real phenopackets with a known single causal gene, from
  **Monarch Phenopacket Store** (<https://github.com/monarch-initiative/phenopacket-store>).
  The HPO profile is a curated case report, so the agent's PubMed search can retrieve the
  *source paper the curator annotated from* and recover the gene without phenotype→gene
  reasoning (RESULTS.md caveat #3).
- **`synthetic`** (pilot 2): patients whose HPO profile is a frequency-weighted draw from the
  **HPOA** aggregate annotations for a disease, generated with
  [`phenotype2phenopacket`](https://github.com/monarch-initiative/phenotype2phenopacket)
  (the pheval ecosystem). No single publication underlies the profile, so that shortcut is
  gone. Disease selection matches the pilot tier — same Phenopacket-Store single-gene disease
  universe, with the 10 pilot diseases always included for a direct A/B; see
  `benchmark/synthesize.py`.
- **`synthetic_hard`** (pilot 3): same generator, but built so Exomiser's phenotype-only
  baseline actually **fails** (pilot 2 showed the synthetic baseline near-ceiling, leaving no
  headroom). Rare/recently-described, sparsely-annotated single-gene diseases; sparse profiles;
  and injected off-target distractor HPO terms (drawn from the HPOA universe minus the
  disease's own terms) that mislead Exomiser's phenotype match while keeping the true signal.
  `benchmark.synthesize --hard`.

Two non-obvious requirements the harness enforces:

1. **Strip the diagnosis before submitting.** Source files contain the answer
   (the `interpretations`/`diagnosis` block, plus the source PMID for `store` and the OMIM
   disease for `synthetic`). Exomiser ignores it (it reads only HPO terms), but the LLM agent
   would *read it and cheat*. The harness feeds OS only `id + subject + phenotypicFeatures +
   metaData`, keeps the ground-truth gene separately, and `prepare.py` refuses to write any
   sanitized file in which the gene symbol, source PMID, or OMIM id survives.
2. **Match genes by stable ID** (Entrez / Ensembl / HGNC), not symbol — symbols have aliases
   and would cause false misses.

## Harness flow

```
synthesize.py  (synthetic only) select tier-matched diseases → p2p create + add-genes (HPOA-sampled
               profile + causal gene), resampling below --min-hpo → data/ppkts_synthetic/
prepare.py     load phenopackets → extract (HPO, true gene id) → strip diagnosis → verify no leakage →
               <sanitized>/ + <ground_truth>.csv
run.py         for each: POST /api/v1/jobs (sanitized phenopacket + prompt, max_iterations≈2) → poll →
               pull exomiser_ranking.json + the rerank output → <results>/<id>/
score.py       true-gene rank in baseline vs rerank → hits@k, McNemar, MRR → report + per-case CSV
```

`--corpus store` writes `sanitized/`, `ground_truth.csv`, `results/`; `--corpus synthetic` writes
the `*_synthetic` siblings, so the two never collide.

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
./download.sh          # Phenopacket Store bundle + HPOA / genes_to_disease / hp.obo

# pilot 1 — real Phenopacket Store cases (source-paper retrievable)
uv run python -m benchmark.prepare --corpus store --n 10
uv run python -m benchmark.run     --corpus store
uv run python -m benchmark.score   --corpus store

# pilot 2 — synthetic HPOA-sampled cases (no source paper)
uv sync --group synth                                       # installs phenotype2phenopacket
uv run python -m benchmark.synthesize --n 20 --seed 0      # → data/ppkts_synthetic/
uv run python -m benchmark.prepare    --corpus synthetic   # → sanitized_synthetic/ + ground_truth_synthetic.csv
uv run python -m benchmark.run        --corpus synthetic
uv run python -m benchmark.score      --corpus synthetic

# pilot 3 — hard cases (rare/recent + sparse + noisy; built to make the baseline fail)
uv run python -m benchmark.synthesize --hard --n 20 --seed 0          # → data/ppkts_synthetic_hard/
uv run python -m benchmark.prepare    --corpus synthetic_hard --min-hpo 3
uv run python -m benchmark.run        --corpus synthetic_hard
uv run python -m benchmark.score      --corpus synthetic_hard
```

`OS_API_KEY` and `.env` are gitignored — never commit them.

## Status

Harness implemented for both corpora. Pilot 1 (10 `store` cases) and pilot 2 (20 `synthetic`
cases) are both run and scored — see RESULTS.md.
