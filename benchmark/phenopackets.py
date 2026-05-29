"""Load GA4GH phenopackets from the Phenopacket Store bundle, extract ground truth,
and sanitize (strip the diagnosis) before they go to OpenScientist.

Structure (validated against release 0.1.26):
  phenotypicFeatures[].type.id                         -> HPO terms (skip excluded)
  interpretations[].diagnosis.disease                  -> {id: OMIM:..., label}
  interpretations[].diagnosis.genomicInterpretations[]
    .variantInterpretation.variationDescriptor.geneContext -> {valueId: HGNC:..., symbol}
"""

from __future__ import annotations

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
    """Unique {symbol: hgnc_id} across all genomic interpretations."""
    genes: dict[str, str] = {}
    for interp in ppkt.get("interpretations", []) or []:
        for gi in (interp.get("diagnosis") or {}).get(
            "genomicInterpretations", []
        ) or []:
            gc = (
                (gi.get("variantInterpretation") or {}).get("variationDescriptor") or {}
            ).get("geneContext") or {}
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


def load_cases(min_hpo: int = 4) -> list[Case]:
    """Eligible cases: exactly one causal gene and >= min_hpo (non-excluded) HPO terms."""
    cases: list[Case] = []
    for f in sorted(config.PPKTS.rglob("*.json")):
        try:
            ppkt = json.loads(f.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        genes = _causal_genes(ppkt)
        hpo = _hpo_ids(ppkt)
        if len(genes) != 1 or len(hpo) < min_hpo:
            continue
        sym, hgnc = next(iter(genes.items()))
        cases.append(
            Case(
                case_id=f"{f.parent.name}__{f.stem}",
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
    """Return the phenopacket with ONLY id/subject/phenotypicFeatures/metaData — the
    interpretations/diseases (the answer) are dropped so the agent can't read it."""
    ppkt = json.loads(ppkt_path.read_text())
    return {
        k: ppkt[k]
        for k in ("id", "subject", "phenotypicFeatures", "metaData")
        if k in ppkt
    }
