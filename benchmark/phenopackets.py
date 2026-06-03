"""Load GA4GH phenopackets from the Phenopacket Store bundle, extract ground truth,
and sanitize (strip the diagnosis) before they go to OpenScientist.

Structure (validated against release 0.1.26):
  phenotypicFeatures[].type.id                         -> HPO terms (skip excluded)
  interpretations[].diagnosis.disease                  -> {id: OMIM:..., label}
  interpretations[].diagnosis.genomicInterpretations[]
    .variantInterpretation.variationDescriptor.geneContext -> {valueId: HGNC:..., symbol}
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

from . import config


@dataclass
class Case:
    case_id: str  # filesystem-safe, unique: "<GENE_cohort>__<file stem>"
    ppkt_id: str
    path: Path
    hpo_ids: list[str]
    gene_symbol: str
    gene_hgnc: str
    omim_disease: str | None


def _hpo_ids(ppkt: dict) -> list[str]:
    out = []
    for f in ppkt.get("phenotypicFeatures", []) or []:
        if f.get("excluded"):
            continue
        hid = (f.get("type") or {}).get("id", "")
        if hid.startswith("HP:"):
            out.append(hid)
    return out


def _causal_genes(ppkt: dict) -> dict[str, str]:
    """Unique {symbol: hgnc_id} across all genomic interpretations.

    Handles both layouts seen in our two corpora: Phenopacket Store nests the gene
    under ``variantInterpretation.variationDescriptor.geneContext``; the synthetic
    phenopackets from ``p2p add-genes`` attach a bare ``gene`` GeneDescriptor.
    """
    genes: dict[str, str] = {}
    for interp in ppkt.get("interpretations", []) or []:
        for gi in (interp.get("diagnosis") or {}).get(
            "genomicInterpretations", []
        ) or []:
            gc = (
                (gi.get("variantInterpretation") or {}).get("variationDescriptor") or {}
            ).get("geneContext") or gi.get("gene") or {}
            sym = gc.get("symbol")
            if sym:
                genes.setdefault(sym, gc.get("valueId") or "")
    return genes


def _omim_disease(ppkt: dict) -> str | None:
    for interp in ppkt.get("interpretations", []) or []:
        d = (interp.get("diagnosis") or {}).get("disease") or {}
        if d.get("id"):
            return d["id"]
    return None


def load_cases(min_hpo: int = 4, ppkts_dir: Path | None = None, flat: bool = False) -> list[Case]:
    """Eligible cases: exactly one causal gene and >= min_hpo (non-excluded) HPO terms.

    ``ppkts_dir`` defaults to the Phenopacket Store bundle. ``flat`` controls the
    case id: the Store nests files under a ``<GENE>/`` dir (id ``<GENE>__<stem>``);
    synthesize writes one flat dir, so we key the id on the gene symbol instead.
    """
    root = ppkts_dir or config.PPKTS
    cases: list[Case] = []
    for f in sorted(root.rglob("*.json")):
        try:
            ppkt = json.loads(f.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        genes = _causal_genes(ppkt)
        hpo = _hpo_ids(ppkt)
        if len(genes) != 1 or len(hpo) < min_hpo:
            continue
        sym, hgnc = next(iter(genes.items()))
        case_id = f"{sym}__{f.stem}" if flat else f"{f.parent.name}__{f.stem}"
        cases.append(
            Case(
                case_id=case_id,
                ppkt_id=ppkt.get("id", f.stem),
                path=f,
                hpo_ids=hpo,
                gene_symbol=sym,
                gene_hgnc=hgnc,
                omim_disease=_omim_disease(ppkt),
            )
        )
    return cases


def sanitize(ppkt_path: Path) -> dict:
    """Return the phenopacket scrubbed of every identifying provenance field so the
    agent cannot look up the answer.

    A 10-case pilot (2026-05-29) hit 100% hits@1 because the agent fetched the
    source paper via PubMed — the leakage was hiding in `id` (carries the source
    PMID), `subject.id` (paper-specific patient label), and `metaData` (especially
    `externalReferences`, which carry the source paper, plus other fields that
    sometimes carry the gene symbol). Only what Exomiser needs to parse and run is
    preserved:

      - a stable but non-identifying `id` (sha256 of the filename, truncated)
      - subject demographics with a generic id (no paper-derived label)
      - phenotypicFeatures (HPO terms — the actual signal)
      - metaData reduced to `phenopacketSchemaVersion` + `resources`
        (ontology version stubs Exomiser needs; everything else dropped)

    `benchmark/prepare.py` re-verifies — refusing to write any sanitized file in
    which the gene symbol or source PMID survives — so silent regression is
    impossible.
    """
    ppkt = json.loads(ppkt_path.read_text())
    case_hash = hashlib.sha256(ppkt_path.name.encode()).hexdigest()[:12]
    out: dict = {"id": f"case_{case_hash}"}

    if "subject" in ppkt:
        sub = {k: v for k, v in ppkt["subject"].items() if k != "id"}
        sub["id"] = f"subject_{case_hash}"
        out["subject"] = sub

    if "phenotypicFeatures" in ppkt:
        out["phenotypicFeatures"] = ppkt["phenotypicFeatures"]

    md = ppkt.get("metaData") or {}
    clean_md: dict = {}
    for k in ("phenopacketSchemaVersion", "resources"):
        if k in md:
            clean_md[k] = md[k]
    out["metaData"] = clean_md

    return out
