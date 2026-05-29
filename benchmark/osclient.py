"""Thin client for the OpenScientist REST API."""

from __future__ import annotations

import time
from pathlib import Path

import httpx

from . import config


def client() -> httpx.Client:
    if not config.OS_API_KEY:
        raise SystemExit("OS_API_KEY is not set (see .env / .env.example).")
    return httpx.Client(
        base_url=config.OS_API_URL,
        headers={"Authorization": f"Bearer {config.OS_API_KEY}"},
        timeout=60.0,
    )


def create_job(
    http: httpx.Client,
    *,
    research_question: str,
    phenopacket_path: Path,
    max_iterations: int,
) -> str:
    """POST /api/v1/jobs with the phenopacket as a data file; return the job id."""
    with phenopacket_path.open("rb") as fh:
        resp = http.post(
            "/api/v1/jobs",
            data={
                "research_question": research_question,
                "max_iterations": str(max_iterations),
            },
            files={"data_files": (phenopacket_path.name, fh, "application/json")},
        )
    resp.raise_for_status()
    return resp.json()["id"]


def get_status(http: httpx.Client, job_id: str) -> dict:
    resp = http.get(f"/api/v1/jobs/{job_id}/status")
    resp.raise_for_status()
    return resp.json()


def wait_for(
    http: httpx.Client, job_id: str, *, poll: float = 10.0, timeout: float = 3600.0
) -> str:
    """Poll until the job reaches a terminal state; return that status (or 'timeout')."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        st = get_status(http, job_id).get("status")
        if st in ("completed", "failed", "cancelled"):
            return st
        time.sleep(poll)
    return "timeout"
