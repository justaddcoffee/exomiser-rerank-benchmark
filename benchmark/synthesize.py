"""Generate synthetic patient phenopackets with phenotype2phenopacket (pheval ecosystem).

Why synthetic: the pilot (RESULTS.md) ran on real Phenopacket Store case reports, so
the agent's PubMed search could retrieve the *source paper the curator annotated from*
and recover the gene without phenotype->gene reasoning. Here each patient's HPO profile
is instead a frequency-weighted draw from the HPOA aggregate annotations for a disease
(`p2p create`), with the causative gene attached from `genes_to_disease.txt`
(`p2p add-genes`). No single publication underlies the profile, so that shortcut is gone.

Disease selection ("match pilot tier"): the candidate pool is the same universe the pilot
drew from -- OMIM diseases present in Phenopacket Store with exactly one causal gene -- kept
only where the fresh HPOA has annotations and `genes_to_disease` maps the disease to exactly
one gene (unambiguous ground truth). The 10 pilot diseases are always included; the rest are
a seeded random sample from the pool, restricted to richly-annotated diseases so the
downsampled profile reliably clears `--min-hpo`.

Reproducibility note: the *disease list* is deterministic (seed). The per-patient HPO
sampling inside p2p uses system randomness, so the exact profiles differ run to run; we
reject and resample any draw below `--min-hpo`. Output goes to the `synthetic` corpus
(`config.SYNTHETIC.ppkts`); feed it to the harness with `prepare/run/score --corpus synthetic`.

Prereqs (./download.sh): data/phenotype.hpoa, data/genes_to_disease.txt, data/hp.obo,
and the Phenopacket Store bundle in data/ppkts/.
"""

from __future__ import annotations

import argparse
import collections
import csv
import json
import random
import subprocess
import sys
import tempfile
from pathlib import Path

from . import config, phenopackets

# The 10 OMIM diseases behind pilot 1 (RESULTS.md); always included so the synthetic
# redo is a direct A/B against the paper-leaky run.
PILOT_OMIMS = [
    "OMIM:151660",  # LMNA
    "OMIM:181350",  # LMNA
    "OMIM:615471",  # FBXL4
    "OMIM:616907",  # CAPN1
    "OMIM:620155",  # SETD2
    "OMIM:620663",  # ERI1
    "OMIM:620849",  # STK33
    "OMIM:620887",  # FDXR
    "OMIM:620988",  # DHX9
    "OMIM:621254",  # ITPR3
]

P2P = Path(sys.executable).with_name("p2p")  # the venv's phenotype2phenopacket CLI


def _hpoa_term_counts() -> collections.Counter:
    """Non-negated phenotype (aspect P) annotation count per OMIM disease."""
    counts: collections.Counter = collections.Counter()
    with config.HPOA.open() as fh:
        for line in fh:
            if line.startswith(("#", "database_id")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("OMIM:") and f[10] == "P" and f[2] != "NOT":
                counts[f[0]] += 1
    return counts


def _g2d_gene_counts() -> dict[str, set[str]]:
    """OMIM disease -> set of gene symbols from genes_to_disease.txt."""
    g2d: dict[str, set[str]] = collections.defaultdict(set)
    with config.GENES_TO_DISEASE.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g2d[row["disease_id"]].add(row["gene_symbol"])
    return g2d


def select_diseases(n: int, seed: int, pool_min_hpo: int, extra_min_hpo: int) -> list[str]:
    """Pilot OMIMs + a seeded sample of tier-matched diseases, total n."""
    store = collections.defaultdict(set)
    for c in phenopackets.load_cases(min_hpo=pool_min_hpo):
        if c.omim_disease:
            store[c.omim_disease].add(c.gene_symbol)
    single_store = {o for o, genes in store.items() if len(genes) == 1}

    hpoa = _hpoa_term_counts()
    g2d = _g2d_gene_counts()

    def eligible(omim: str, hpoa_floor: int) -> bool:
        return (
            omim in single_store
            and hpoa.get(omim, 0) >= hpoa_floor
            and len(g2d.get(omim, set())) == 1
        )

    chosen = [o for o in PILOT_OMIMS if eligible(o, pool_min_hpo)]
    missing = [o for o in PILOT_OMIMS if o not in chosen]
    if missing:
        print(f"warning: pilot OMIMs not eligible in current data, skipped: {missing}")

    # Backfill from richly-annotated tier-matched diseases (robust to downsampling).
    extra_pool = sorted(
        o for o in single_store if eligible(o, extra_min_hpo) and o not in chosen
    )
    rng = random.Random(seed)
    while len(chosen) < n and extra_pool:
        chosen.append(extra_pool.pop(rng.randrange(len(extra_pool))))
    if len(chosen) < n:
        print(f"warning: pool exhausted; selected {len(chosen)} of requested {n}")
    return sorted(set(chosen))


def _p2p(*args: str) -> None:
    subprocess.run([str(P2P), *args], check=True)


def _non_excluded_hpo(ppkt_path: Path) -> int:
    ppkt = json.loads(ppkt_path.read_text())
    return len(phenopackets._hpo_ids(ppkt))


def generate(omims: list[str], min_hpo: int, max_rounds: int) -> dict[str, Path]:
    """create + add-genes per disease, resampling any profile below min_hpo.

    Returns {omim: best phenopacket path} written into the synthetic corpus dir.
    """
    out_dir = config.SYNTHETIC.ppkts
    out_dir.mkdir(parents=True, exist_ok=True)
    best: dict[str, tuple[int, Path]] = {}  # omim -> (n_hpo, path)
    remaining = list(omims)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_root = Path(tmp)
        for rnd in range(1, max_rounds + 1):
            if not remaining:
                break
            print(f"[round {rnd}] sampling {len(remaining)} disease(s)...", flush=True)
            ids_file = tmp_root / f"ids_{rnd}.txt"
            ids_file.write_text("\n".join(remaining) + "\n")
            created = tmp_root / f"create_{rnd}"
            withgenes = tmp_root / f"genes_{rnd}"
            cache = ["-c", str(config.HP_OBO)] if config.HP_OBO.exists() else []
            _p2p("create", "-p", str(config.HPOA), "-l", str(ids_file),
                 "-o", str(created), *cache)
            _p2p("add-genes", "-p", str(created), "-g", str(config.GENES_TO_DISEASE),
                 "-o", str(withgenes), "-i", "hgnc_id")

            still: list[str] = []
            for omim in remaining:
                src = withgenes / f"{omim.replace(':', '_')}_patient_1.json"
                if not src.exists():
                    still.append(omim)
                    continue
                n = _non_excluded_hpo(src)
                dst = out_dir / src.name
                if n >= min_hpo:
                    dst.write_text(src.read_text())
                    best[omim] = (n, dst)
                else:
                    if omim not in best or n > best[omim][0]:
                        dst.write_text(src.read_text())
                        best[omim] = (n, dst)
                    still.append(omim)
            remaining = still

    weak = {o: n for o, (n, _) in best.items() if n < min_hpo}
    if weak:
        print(f"warning: below --min-hpo after {max_rounds} rounds (kept best): {weak}")
    return {o: p for o, (_, p) in best.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=20, help="number of synthetic cases")
    ap.add_argument("--seed", type=int, default=0, help="disease-selection seed")
    ap.add_argument(
        "--min-hpo", type=int, default=4, help="minimum non-excluded HPO terms per patient"
    )
    ap.add_argument(
        "--pool-min-hpo",
        type=int,
        default=4,
        help="HPOA annotation floor for the candidate pool (incl. pilot diseases)",
    )
    ap.add_argument(
        "--extra-min-hpo",
        type=int,
        default=15,
        help="HPOA annotation floor for backfill diseases (robust to downsampling)",
    )
    ap.add_argument("--max-rounds", type=int, default=12, help="resampling rounds")
    args = ap.parse_args()

    for p in (config.HPOA, config.GENES_TO_DISEASE):
        if not p.exists():
            raise SystemExit(f"{p} not found — run ./download.sh first.")
    if not config.PPKTS.exists():
        raise SystemExit(
            f"{config.PPKTS} not found — run ./download.sh (Phenopacket Store) first."
        )

    omims = select_diseases(args.n, args.seed, args.pool_min_hpo, args.extra_min_hpo)
    print(f"selected {len(omims)} diseases:\n  " + "\n  ".join(omims))

    config.SYNTHETIC.ppkts.mkdir(parents=True, exist_ok=True)
    for stale in config.SYNTHETIC.ppkts.glob("*.json"):
        stale.unlink()

    written = generate(omims, args.min_hpo, args.max_rounds)
    print(
        f"\nwrote {len(written)} synthetic phenopackets -> {config.SYNTHETIC.ppkts}/\n"
        f"next: uv run python -m benchmark.prepare --corpus synthetic"
    )


if __name__ == "__main__":
    main()
