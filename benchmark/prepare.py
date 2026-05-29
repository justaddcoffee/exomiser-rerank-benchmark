"""Build the case set: fetch phenopackets, extract ground truth, strip the diagnosis.

Source: Monarch Phenopacket Store (github.com/monarch-initiative/phenopacket-store) — real
GA4GH phenopackets with a curated diagnosis (causal gene + OMIM disease).

For each selected case this produces:
  - sanitized/<id>.json : phenopacket with ONLY id + subject + phenotypicFeatures + metaData
                          (the interpretations/diagnosis block is removed so the LLM can't
                          read the answer; Exomiser only needs the HPO terms anyway).
  - ground_truth.csv    : id, gene_symbol, gene_entrez, gene_ensembl, gene_hgnc, omim_disease

Selection favors single-gene Mendelian cases; the pilot takes N=10 (diverse genes).

TODO(impl):
  - fetch_phenopackets(n): pull N phenopackets (raw JSON from the store).
  - extract_truth(ppkt): pull (HPO ids, causal gene with stable IDs, OMIM disease) from
    interpretations -> diagnosis -> genomicInterpretations -> variationDescriptor/geneContext.
  - sanitize(ppkt): drop interpretations/diagnosis; keep id/subject/phenotypicFeatures/metaData.
"""

from __future__ import annotations

import argparse


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, default=10, help="number of phenopackets (pilot=10)")
    parser.add_argument("--out", default="sanitized", help="output dir for stripped phenopackets")
    args = parser.parse_args()
    raise NotImplementedError(
        f"prepare: fetch {args.n} phenopackets -> strip diagnosis -> {args.out}/ + ground_truth.csv"
    )


if __name__ == "__main__":
    main()
