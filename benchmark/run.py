"""Submit each sanitized phenopacket to OpenScientist and collect the two rankings.

Client of the OS REST API (OS_API_URL / OS_API_KEY from .env). Per case:
  1. POST /api/v1/jobs  (multipart: research_question=PROMPT, data_files=<sanitized phenopacket>,
     max_iterations=2)
  2. poll GET /api/v1/jobs/{id}/status until completed/failed
  3. retrieve the job's `exomiser_ranking.json` (baseline) and the reranked gene list, save to
     results/<id>/.

Concurrency capped at OS_MAX_INFLIGHT. The prompt below defines the rerank output schema
explicitly, since OS treats reranking as opt-in and does not mandate a format.

TODO(impl):
  - submit(client, ppkt_path) -> job_id
  - poll(client, job_id) -> status
  - collect(job_id) -> {"exomiser_ranking": ..., "reranked": ...}  (how to fetch artifacts —
    API endpoint vs job-dir on disk — depends on where OS runs; local pilot can read the job dir).
"""

from __future__ import annotations

PROMPT = """\
Attached is a patient phenopacket containing HPO phenotype terms (no diagnosis).

1. Run Exomiser in phenotype-only mode (the `run_exomiser` tool, preset=phenotype-only) to
   produce a ranked differential of candidate GENES. This emits `exomiser_ranking.json` — the
   immutable baseline; do not modify it.
2. Re-rank the top candidate genes using the evidence resources available to you (PubMed,
   ClinVar, OMIM, the database skills). Every reorder must cite at least one source.
3. Write the re-ranking to `reranked.json` as a JSON array of objects, best first:
   {"rank": <int>, "geneSymbol": <str>, "geneId": <stable id e.g. ENSG/Entrez/HGNC>,
    "exomiserRank": <int>, "rationale": <str>,
    "evidence": [{"source": <str>, "id": <str>, "note": <str>}]}
4. Report the single top candidate gene.
"""


def main() -> None:
    raise NotImplementedError("run: submit sanitized phenopackets via the OS API + collect rankings")


if __name__ == "__main__":
    main()
