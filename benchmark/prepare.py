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

from . import config, phenopackets


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=10, help="number of cases (pilot=10)")
    ap.add_argument(
        "--min-hpo", type=int, default=4, help="minimum non-excluded HPO terms"
    )
    ap.add_argument("--seed", type=int, default=0, help="sampling seed (reproducible)")
    args = ap.parse_args()

    if not config.PPKTS.exists():
        raise SystemExit(
            f"{config.PPKTS} not found — download the Phenopacket Store bundle first."
        )

    cases = phenopackets.load_cases(min_hpo=args.min_hpo)
    if not cases:
        raise SystemExit("no eligible cases found")
    random.Random(args.seed).shuffle(cases)
    chosen = cases[: args.n]

    config.SANITIZED.mkdir(parents=True, exist_ok=True)
    rows = []
    for c in chosen:
        (config.SANITIZED / f"{c.case_id}.json").write_text(
            json.dumps(phenopackets.sanitize(c.path), indent=2)
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

    with config.GROUND_TRUTH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    print(
        f"prepared {len(rows)} cases (of {len(cases)} eligible) "
        f"-> {config.SANITIZED}/ + {config.GROUND_TRUTH.name}"
    )


if __name__ == "__main__":
    main()
