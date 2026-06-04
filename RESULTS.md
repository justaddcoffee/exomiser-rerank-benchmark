# Pilot 1 — 10 Phenopacket-Store cases (2026-05)

A first proof-of-concept run of [OpenScientist](https://github.com/openscientist-io/openscientist) with the [Exomiser](https://www.sanger.ac.uk/tool/exomiser/) MCP tool, asking: **does an LLM-driven evidence rerank of Exomiser's phenotype-only candidate list improve gene prioritization on real case-report phenopackets?**

## Setup

10 single-causal-gene cases (≥4 non-excluded HPO terms; reproducible seed=0) drawn from the GA4GH [Phenopacket Store](https://github.com/monarch-initiative/phenopacket-store). Each phenopacket is sanitized before it reaches the agent — the diagnosis (`interpretations`, `diseases`), the source paper `id`/PMID, the gene symbol in `metaData`, and any other identifying provenance are stripped; `benchmark/prepare.py` refuses to write any sanitized file in which the true gene symbol or the source PMID survives, so the agent receives HPO terms + demographics only.

Each case becomes one OpenScientist job with this prompt:

> 1. Run Exomiser in **phenotype-only mode** (`run_exomiser` tool) to produce `exomiser_ranking.json` — the immutable baseline.
> 2. **Re-rank the top candidates** using cited evidence (PubMed, ClinVar, OMIM, OS's database skills). Every reorder must cite at least one source.
> 3. Write the rerank to `reranked.json` as an array of `{rank, geneSymbol, geneId, exomiserRank, rationale, evidence[]}`.

2 LLM iterations per job; 10 jobs run sequentially on a single Mac (emulated amd64 Exomiser via the agent container). The harness then collects both rankings, computes the true causal gene's position by symbol match in each, and reports hits@{1, 3, 5, 10}, MRR, and McNemar's exact paired test.

## Result

| k  | Exomiser baseline (phenotype-only) | OS-reranked | McNemar p |
|----|---|---|---|
| 1  | 3/10 (30 %)  | **10/10 (100 %)** | **0.016** |
| 3  | 4/10 (40 %)  | 10/10 (100 %)     | 0.031 |
| 5  | 6/10 (60 %)  | 10/10 (100 %)     | 0.125 |
| 10 | 6/10 (60 %)  | 10/10 (100 %)     | 0.125 |

MRR: baseline **0.41** → reranked **1.00**. Reranking moved the true gene up in 7/10 cases, down in 0.

Cases: `LMNA__PMID_37843397_TR36`, `ITPR3__PMID_39560673_P2`, `DHX9__PMID_37467750_Individual_12`, `FDXR__PMID_29040572_Family_9_individual`, `SETD2__PMID_37372360_P18`, `ERI1__PMID_37352860_Individual_1A`, `FBXL4__PMID_28940506_Family_28_Individual_34`, `LMNA__PMID_10939567_Patient_Spo11`, `CAPN1__PMID_28321562_sister`, `STK33__PMID_34155512_III_1`.

## Caveats (important for reading the number)

1. **n = 10.** The McNemar p-values are a proof-of-concept signal, not a study finding.
2. **"Famous monogenic" tier.** These genes have rich literature with characteristic phenotype-gene associations that an LLM with PubMed is well-equipped to recover. Rare or recently-discovered variants would be a substantially tougher test and the next thing to add to the set.
3. **Source-paper retrieval, not source-paper leakage.** The input is now clean of the source PMID after sanitization — but in 9/10 cases the agent's `search_pubmed` still organically retrieves the source case-report paper, because for these phenotypes that paper *is* the canonical reference. That's not leakage from the input — but it does mean some of the 100 % hits@1 is "the agent found the same paper the curator did" rather than purely "the agent reasoned from HPO terms to gene." A cleaner next experiment would block the per-case source PMID from `search_pubmed`'s retrieval set (small change in OS).

## Reproducing

```bash
git clone https://github.com/justaddcoffee/exomiser-rerank-benchmark.git
cd exomiser-rerank-benchmark
./download.sh                                 # pulls Phenopacket Store release 0.1.26
cp .env.example .env && $EDITOR .env          # OS_API_URL, OS_API_KEY, OS_JOBS_DIR
uv run python -m benchmark.prepare --n 10     # selects 10 cases, writes sanitized/ + ground_truth.csv
uv run python -m benchmark.run                # submits each to OS, polls, collects rankings
uv run python -m benchmark.score              # hits@k, McNemar, MRR
```

Output:
- `results/<case_id>/{exomiser_ranking.json,reranked.json,meta.json}` — per-case artifacts collected from the OS job dir.
- `results/scores.csv` — per-case true-gene ranks in baseline vs reranked.

# Pilot 2 — 20 synthetic cases (2026-06)

The big asterisk on pilot 1 was caveat #3: with real case-report phenopackets, the agent's
`search_pubmed` could retrieve *the very paper the curator annotated from* and recover the gene
without reasoning from phenotype to gene. Pilot 2 removes that shortcut at the root.

## Setup

20 synthetic patients generated with [`phenotype2phenopacket`](https://github.com/monarch-initiative/phenotype2phenopacket)
(the pheval ecosystem; `benchmark/synthesize.py`). Each patient's HPO profile is a
frequency-weighted draw from the **HPOA aggregate** annotations for an OMIM disease — **no single
publication underlies the profile**, so the source-paper shortcut is gone. The causal gene is
attached from HPO's `genes_to_disease`. Disease selection matches the pilot tier (same
Phenopacket-Store single-gene disease universe); the 10 pilot-1 diseases are always included, so
10 of the 20 cases are a direct A/B against pilot 1. Sanitization is unchanged, and `prepare.py`
now also refuses any sanitized file in which the OMIM id survives. Same prompt, same 2 LLM
iterations, same local emulated Exomiser (phenotype-only, data `2512_phenotype`). One job hit the
per-job poll timeout but still returned both rankings; 0 failures.

## Result

| k  | Exomiser baseline (phenotype-only) | OS-reranked | McNemar p |
|----|---|---|---|
| 1  | 15/20 (75 %) | **20/20 (100 %)** | 0.062 |
| 3  | 17/20 (85 %) | 20/20 (100 %)     | 0.250 |
| 5  | 17/20 (85 %) | 20/20 (100 %)     | 0.250 |
| 10 | 19/20 (95 %) | 20/20 (100 %)     | 1.000 |

MRR: baseline **0.81** → reranked **1.00**. Reranking moved the true gene up in 5/20 cases, down
in 0. The five with headroom — the only baseline misses — were all recovered: GNAS 11→1,
ITPR3 10→1, STK33 10→1, TUBA4A 3→1, IRF4 2→1.

**Headline: the reranking gain survives removal of the source-paper shortcut.** With no paper to
retrieve, OS still reaches 100 % hits@1 (20/20). On the 10 shared diseases, reranked rank = 1 in
*both* pilots; e.g. CAPN1 (P1 5→1 / P2 1→1), FBXL4 (16→1 / 1→1), SETD2 (25→1 / 1→1), STK33
(2→1 / 10→1).

## The caveat, now inverted

The thing that changed between pilots is the **baseline**, not the rerank. On the 10 shared
diseases, baseline hits@1 jumped **3/10 → 8/10**; overall pilot-2 baseline is 75 % @1 vs pilot-1's
30 %. Synthetic HPOA-aggregate profiles are more *prototypical* than messy real case reports, so
Exomiser's phenotype-only match is much stronger and there is little headroom left. That — not a
weaker reranker — is why the McNemar p-values aren't significant here: only 5 cases could move,
and all 5 did, with 0 harm. (STK33 actually got *harder* synthetically, baseline 2→10, and
reranking still recovered it.)

So pilot 2 cleanly retires the leakage worry, and shifts the open question: to stress-test
reranking you now need cases where the baseline genuinely **fails** — rarer/recently-described
diseases, or noisier/sparser profiles (lower `--min-hpo`, more added noise).

## Reproducing

```bash
./download.sh                                          # + HPOA / genes_to_disease / hp.obo
uv sync --group synth                                  # phenotype2phenopacket
uv run python -m benchmark.synthesize --n 20 --seed 0  # disease list is deterministic;
                                                       # per-patient HPO sampling is stochastic
uv run python -m benchmark.prepare --corpus synthetic  # sanitized_synthetic/ + ground_truth_synthetic.csv
uv run python -m benchmark.run     --corpus synthetic
uv run python -m benchmark.score   --corpus synthetic  # -> results_synthetic/scores.csv
```

## What's next

- **~100-case run** with source-PMID retrieval blocked, balanced across difficulty tiers.
- Per-case **headroom check** — for the cases where the baseline already nails it, does reranking ever *hurt*? (Pilot says no — 0/10 went down — but n is small.)
- **Evidence quality** as a separate axis: even when reranking lands the right gene, is the cited evidence load-bearing for the decision, or does the rationale stand without it? Hook for [openscientist-io/openscientist#191](https://github.com/openscientist-io/openscientist/pull/191) (the citation-validation work).

— j. reese, 2026-06
