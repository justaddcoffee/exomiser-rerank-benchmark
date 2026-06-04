"""Submit each sanitized phenopacket to OpenScientist and collect both rankings.

Per case: POST a discovery job (phenopacket + PROMPT, phenotype-only Exomiser + rerank),
wait for completion, then copy `exomiser_ranking.json` (baseline) and `reranked.json`
(agent rerank) out of the OS job dir into results/<case_id>/.

Collection reads OS_JOBS_DIR (local OS instance). Jobs are processed sequentially — fine
for the 10-case pilot; parallelize up to OS_MAX_INFLIGHT for the full run.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil

from . import config, osclient

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


def _collect(job_id: str, case_id: str, cp: config.Corpus) -> dict:
    """Copy the baseline + rerank rankings from the OS job dir into <results>/<case_id>/."""
    out = cp.results / case_id
    out.mkdir(parents=True, exist_ok=True)
    found: dict[str, bool] = {}
    if config.OS_JOBS_DIR is None:
        out.joinpath("meta.json").write_text(
            json.dumps({"job_id": job_id, "error": "OS_JOBS_DIR unset"})
        )
        return found
    jd = config.OS_JOBS_DIR / job_id
    # exomiser_ranking.json lives under data/exomiser_results/<hex>/; reranked.json wherever
    # the agent wrote it — search broadly to be robust.
    for name in ("exomiser_ranking.json", "reranked.json"):
        hits = sorted(jd.rglob(name))
        if hits:
            shutil.copy(hits[-1], out / name)
            found[name] = True
    out.joinpath("meta.json").write_text(
        json.dumps({"job_id": job_id, "found": found}, indent=2)
    )
    return found


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--corpus",
        choices=sorted(config.CORPORA),
        default="store",
        help="which prepared corpus to submit",
    )
    ap.add_argument("--max-iterations", type=int, default=2)
    args = ap.parse_args()

    cp = config.corpus(args.corpus)
    if not cp.ground_truth.exists():
        raise SystemExit(
            f"{cp.ground_truth.name} not found — run "
            f"`python -m benchmark.prepare --corpus {cp.name}` first."
        )
    cases = list(csv.DictReader(cp.ground_truth.open()))
    http = osclient.client()

    for i, c in enumerate(cases, 1):
        cid = c["case_id"]
        ppkt = cp.sanitized / f"{cid}.json"
        job_id = osclient.create_job(
            http,
            research_question=PROMPT,
            phenopacket_path=ppkt,
            max_iterations=args.max_iterations,
        )
        print(f"[{cp.name} {i}/{len(cases)}] {cid}: submitted job {job_id}", flush=True)
        status = osclient.wait_for(http, job_id)
        found = _collect(job_id, cid, cp)
        print(f"            -> {status}; collected {sorted(found)}", flush=True)


if __name__ == "__main__":
    main()
