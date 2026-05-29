"""Paths + environment config (loaded from .env)."""

from __future__ import annotations

import os
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

# OpenScientist REST API
OS_API_URL = os.environ.get("OS_API_URL", "http://localhost:8080").rstrip("/")
OS_API_KEY = os.environ.get("OS_API_KEY", "")

# Where OS writes job artifacts on disk (HOST_PROJECT_DIR/jobs). Needed to collect
# exomiser_ranking.json / reranked.json for a LOCAL OS instance. (Remote OS would
# need an artifact API or SSH instead.)
OS_JOBS_DIR = Path(os.environ["OS_JOBS_DIR"]) if os.environ.get("OS_JOBS_DIR") else None

OS_MAX_INFLIGHT = int(os.environ.get("OS_MAX_INFLIGHT", "2"))
