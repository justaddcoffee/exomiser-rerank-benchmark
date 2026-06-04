"""Paths + environment config (loaded from .env)."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

from dotenv import load_dotenv

load_dotenv()

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
PPKTS = (
    DATA / "ppkts"
)  # unzipped Phenopacket Store bundle (data/ppkts/<release>/<GENE>/*.json)
SANITIZED = REPO / "sanitized"  # diagnosis-stripped phenopackets submitted to OS
RESULTS = REPO / "results"  # collected rankings + scores
GROUND_TRUTH = REPO / "ground_truth.csv"

# Inputs for synthetic-case generation (phenotype2phenopacket); fetched by download.sh.
HPOA = DATA / "phenotype.hpoa"
GENES_TO_DISEASE = DATA / "genes_to_disease.txt"
HP_OBO = DATA / "hp.obo"


@dataclass(frozen=True)
class Corpus:
    """A set of cases with its own input/output locations.

    Two corpora share one harness (prepare/run/score) via the ``--corpus`` flag:
      - ``store``: real GA4GH Phenopacket Store case reports (pilot 1). The HPO
        profile is the curated case, so the source paper can leak via PubMed.
      - ``synthetic``: patients whose HPO profile is sampled from the HPOA
        aggregate by ``benchmark.synthesize`` — no single source paper exists.
    ``flat`` distinguishes the Phenopacket Store's ``<GENE>/<file>`` layout (case
    id = ``<GENE>__<stem>``) from synthesize's single output dir (id =
    ``<gene>__<stem>``, gene resolved from the phenopacket).
    """

    name: str
    ppkts: Path
    sanitized: Path
    ground_truth: Path
    results: Path
    flat: bool


STORE = Corpus("store", PPKTS, SANITIZED, GROUND_TRUTH, RESULTS, flat=False)
SYNTHETIC = Corpus(
    "synthetic",
    DATA / "ppkts_synthetic",
    REPO / "sanitized_synthetic",
    REPO / "ground_truth_synthetic.csv",
    REPO / "results_synthetic",
    flat=True,
)
# Like `synthetic`, but deliberately hard: rare/recent sparsely-annotated diseases,
# sparse profiles, and injected distractor HPO terms -- built so Exomiser's
# phenotype-only baseline actually fails, giving reranking measurable headroom.
SYNTHETIC_HARD = Corpus(
    "synthetic_hard",
    DATA / "ppkts_synthetic_hard",
    REPO / "sanitized_synthetic_hard",
    REPO / "ground_truth_synthetic_hard.csv",
    REPO / "results_synthetic_hard",
    flat=True,
)
CORPORA = {c.name: c for c in (STORE, SYNTHETIC, SYNTHETIC_HARD)}


def corpus(name: str) -> Corpus:
    try:
        return CORPORA[name]
    except KeyError:
        raise SystemExit(
            f"unknown corpus '{name}'; choose from {sorted(CORPORA)}"
        ) from None

# OpenScientist REST API
OS_API_URL = os.environ.get("OS_API_URL", "http://localhost:8080").rstrip("/")
OS_API_KEY = os.environ.get("OS_API_KEY", "")

# Where OS writes job artifacts on disk (HOST_PROJECT_DIR/jobs). Needed to collect
# exomiser_ranking.json / reranked.json for a LOCAL OS instance. (Remote OS would
# need an artifact API or SSH instead.)
OS_JOBS_DIR = Path(os.environ["OS_JOBS_DIR"]) if os.environ.get("OS_JOBS_DIR") else None

OS_MAX_INFLIGHT = int(os.environ.get("OS_MAX_INFLIGHT", "2"))
