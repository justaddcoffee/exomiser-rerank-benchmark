"""Select cases from the Phenopacket Store bundle, strip diagnoses, write ground truth.

Produces:
  sanitized/<case_id>.json  -- phenopacket with the diagnosis removed (fed to OS)
  ground_truth.csv          -- case_id, ppkt_id, gene_symbol, gene_hgnc, omim_disease, n_hpo

Run `download.sh` (or fetch the release zip into data/ppkts/) first.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
import re

from . import config, phenopackets


def _check_no_leakage(sanitized: dict, case: phenopackets.Case) -> None:
    """Refuse to write a sanitized phenopacket in which the answer survives.

    Catches both the obvious case (gene symbol appearing anywhere) and the subtle
    one that broke the first pilot (the source PMID present so the agent can fetch
    the paper). Word-boundary matching on the gene symbol avoids false positives
    on substrings inside HPO labels.
    """
    text = json.dumps(sanitized)
    leaks: list[str] = []
    if case.gene_symbol and re.search(rf"\b{re.escape(case.gene_symbol)}\b", text):
        leaks.append(f"gene_symbol='{case.gene_symbol}'")
    if case.ppkt_id.startswith("PMID_"):
        parts = case.ppkt_id.split("_")
        if len(parts) >= 2 and parts[1].isdigit() and parts[1] in text:
            leaks.append(f"PMID='{parts[1]}'")
    # Synthetic cases carry the diagnosis as an OMIM id (no PMID); make sure neither
    # the id nor its bare number survives so the agent can't shortcut to the disease.
    if case.omim_disease and (
        case.omim_disease in text or case.omim_disease.split(":")[-1] in text
    ):
        leaks.append(f"omim='{case.omim_disease}'")
    if leaks:
        raise SystemExit(
            f"LEAKAGE in sanitized phenopacket for case '{case.case_id}': "
            f"{', '.join(leaks)}. Strengthen phenopackets.sanitize() before running."
        )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--corpus",
        choices=sorted(config.CORPORA),
        default="store",
        help="which corpus to prepare (store=Phenopacket Store, synthetic=HPOA-sampled)",
    )
    ap.add_argument(
        "--n",
        type=int,
        default=0,
        help="number of cases (0 = all eligible; pilot used 10)",
    )
    ap.add_argument(
        "--min-hpo", type=int, default=4, help="minimum non-excluded HPO terms"
    )
    ap.add_argument("--seed", type=int, default=0, help="sampling seed (reproducible)")
    args = ap.parse_args()

    cp = config.corpus(args.corpus)
    if not cp.ppkts.exists():
        hint = (
            "download the Phenopacket Store bundle first (./download.sh)"
            if cp is config.STORE
            else "generate it first (python -m benchmark.synthesize)"
        )
        raise SystemExit(f"{cp.ppkts} not found — {hint}.")

    cases = phenopackets.load_cases(
        min_hpo=args.min_hpo, ppkts_dir=cp.ppkts, flat=cp.flat
    )
    if not cases:
        raise SystemExit("no eligible cases found")
    if args.n:
        random.Random(args.seed).shuffle(cases)
        chosen = cases[: args.n]
    else:
        chosen = cases

    cp.sanitized.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in chosen:
        sanitized = phenopackets.sanitize(c.path)
        _check_no_leakage(sanitized, c)
        (cp.sanitized / f"{c.case_id}.json").write_text(
            json.dumps(sanitized, indent=2)
        )
        rows.append(
            {
                "case_id": c.case_id,
                "ppkt_id": c.ppkt_id,
                "gene_symbol": c.gene_symbol,
                "gene_hgnc": c.gene_hgnc,
                "omim_disease": c.omim_disease or "",
                "n_hpo": len(c.hpo_ids),
            }
        )

    with cp.ground_truth.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"[{cp.name}] prepared {len(rows)} cases (of {len(cases)} eligible) "
        f"-> {cp.sanitized}/ + {cp.ground_truth.name}"
    )


if __name__ == "__main__":
    main()
