"""Generate synthetic patient phenopackets with phenotype2phenopacket (pheval ecosystem).

Why synthetic: the pilot (RESULTS.md) ran on real Phenopacket Store case reports, so
the agent's PubMed search could retrieve the *source paper the curator annotated from*
and recover the gene without phenotype->gene reasoning. Here each patient's HPO profile
is instead a frequency-weighted draw from the HPOA aggregate annotations for a disease
(`p2p create`), with the causative gene attached from `genes_to_disease.txt`
(`p2p add-genes`). No single publication underlies the profile, so that shortcut is gone.

Two presets, selected by `--corpus` (or `--hard`):

  synthetic       -- tier-matched to pilot 1 (same Phenopacket-Store single-gene disease
                     universe; the 10 pilot diseases always included for a direct A/B).
  synthetic_hard  -- deliberately hard, built so Exomiser's phenotype-only baseline *fails*
                     (pilot 2 showed the synthetic baseline is near-ceiling, leaving no
                     headroom): rare/recent + sparsely-annotated diseases, sparse profiles
                     (`--min/--max-hpo`), and injected off-target distractor terms
                     (`--noise-terms`) drawn from the HPOA universe minus the disease's own
                     terms. All true terms are kept, so an evidence-reasoning agent can still
                     get there; Exomiser's raw phenotype match is degraded.

Disease selection ("match pilot tier"): the candidate pool is the same universe the pilot
drew from -- OMIM diseases present in Phenopacket Store with exactly one causal gene -- kept
only where the fresh HPOA has annotations and `genes_to_disease` maps the disease to exactly
one gene (unambiguous ground truth). `--rare` narrows to recently-described (`--recent-omim`),
sparsely-annotated (`--hpoa-min/max`) diseases and drops the well-characterized pilot set.

Reproducibility note: the *disease list* is deterministic (seed). The per-patient HPO
sampling inside p2p uses system randomness, so the exact profiles differ run to run; we
reject and resample any draw outside the `[--min-hpo, --max-hpo]` window. Injected noise is
seeded. Output goes to the chosen corpus (`config.<CORPUS>.ppkts`); feed it to the harness
with `prepare/run/score --corpus <name>`.

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

# The 10 OMIM diseases behind pilot 1 (RESULTS.md); always included by the default preset
# so the synthetic redo is a direct A/B against the paper-leaky run.
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


def _omim_num(omim: str) -> int:
    try:
        return int(omim.split(":")[1])
    except (IndexError, ValueError):
        return 0


def _hpoa_disease_terms() -> dict[str, set[str]]:
    """OMIM disease -> set of its non-negated phenotype (aspect P) HPO term ids."""
    terms: dict[str, set[str]] = collections.defaultdict(set)
    with config.HPOA.open() as fh:
        for line in fh:
            if line.startswith(("#", "database_id")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0].startswith("OMIM:") and f[10] == "P" and f[2] != "NOT":
                terms[f[0]].add(f[3])
    return terms


def _hp_labels() -> dict[str, str]:
    """HP term id -> primary name, parsed from data/hp.obo."""
    labels: dict[str, str] = {}
    cur: str | None = None
    with config.HP_OBO.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "[Term]":
                cur = None
            elif line.startswith("id: HP:"):
                cur = line[4:].strip()
            elif line.startswith("name:") and cur:
                labels[cur] = line[5:].strip()
                cur = None
    return labels


def _g2d_gene_counts() -> dict[str, set[str]]:
    """OMIM disease -> set of gene symbols from genes_to_disease.txt."""
    g2d: dict[str, set[str]] = collections.defaultdict(set)
    with config.GENES_TO_DISEASE.open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            g2d[row["disease_id"]].add(row["gene_symbol"])
    return g2d


def select_diseases(
    n: int,
    seed: int,
    pool_min_hpo: int,
    extra_min_hpo: int,
    *,
    rare: bool = False,
    hpoa_min: int = 6,
    hpoa_max: int = 15,
    recent_omim: int = 619000,
) -> list[str]:
    """Pick n OMIM diseases. Default: pilot OMIMs + a seeded tier-matched sample.
    ``rare``: recently-described (OMIM >= recent_omim), sparsely-annotated
    (hpoa_min..hpoa_max) single-gene diseases, with the pilot set excluded."""
    disease_terms = _hpoa_disease_terms()
    hpoa = {o: len(t) for o, t in disease_terms.items()}
    g2d = _g2d_gene_counts()

    store = collections.defaultdict(set)
    for c in phenopackets.load_cases(min_hpo=1):
        if c.omim_disease:
            store[c.omim_disease].add(c.gene_symbol)
    single_store = {o for o, genes in store.items() if len(genes) == 1}

    def single_gene(omim: str) -> bool:
        return omim in single_store and len(g2d.get(omim, set())) == 1

    rng = random.Random(seed)
    if rare:
        pool = sorted(
            o
            for o in single_store
            if single_gene(o)
            and hpoa_min <= hpoa.get(o, 0) <= hpoa_max
            and _omim_num(o) >= recent_omim
            and o not in PILOT_OMIMS
        )
        if len(pool) < n:
            raise SystemExit(
                f"rare pool has only {len(pool)} diseases (<{n}); widen "
                f"--hpoa-min/max or lower --recent-omim."
            )
        chosen = [pool.pop(rng.randrange(len(pool))) for _ in range(n)]
        return sorted(chosen)

    chosen = [o for o in PILOT_OMIMS if single_gene(o) and hpoa.get(o, 0) >= pool_min_hpo]
    missing = [o for o in PILOT_OMIMS if o not in chosen]
    if missing:
        print(f"warning: pilot OMIMs not eligible in current data, skipped: {missing}")
    extra_pool = sorted(
        o
        for o in single_store
        if single_gene(o) and hpoa.get(o, 0) >= extra_min_hpo and o not in chosen
    )
    while len(chosen) < n and extra_pool:
        chosen.append(extra_pool.pop(rng.randrange(len(extra_pool))))
    if len(chosen) < n:
        print(f"warning: pool exhausted; selected {len(chosen)} of requested {n}")
    return sorted(set(chosen))


def _p2p(*args: str) -> None:
    subprocess.run([str(P2P), *args], check=True)


def _add_noise(path: Path, omim: str, n_noise: int, seed: int) -> int:
    """Append n_noise off-target distractor HPO terms to a phenopacket in place.

    Distractors are drawn from the HPOA phenotype universe minus the target disease's
    own annotated terms (and the terms already present), so they actively mislead
    Exomiser's phenotype match while leaving the true signal intact. Returns the number
    actually added (labels resolved from hp.obo). Module-level caches keep it cheap.
    """
    if n_noise <= 0:
        return 0
    universe = _add_noise.universe  # type: ignore[attr-defined]
    labels = _add_noise.labels  # type: ignore[attr-defined]
    disease_terms = _add_noise.disease_terms  # type: ignore[attr-defined]

    ppkt = json.loads(path.read_text())
    present = {
        (f.get("type") or {}).get("id")
        for f in ppkt.get("phenotypicFeatures", []) or []
    }
    exclude = present | disease_terms.get(omim, set())
    candidates = [t for t in universe if t not in exclude and t in labels]
    rng = random.Random(f"{seed}:{omim}")
    rng.shuffle(candidates)
    added = candidates[:n_noise]
    ppkt.setdefault("phenotypicFeatures", []).extend(
        {"type": {"id": t, "label": labels[t]}} for t in added
    )
    path.write_text(json.dumps(ppkt, indent=2))
    return len(added)


def _non_excluded_hpo(ppkt_path: Path) -> int:
    ppkt = json.loads(ppkt_path.read_text())
    return len(phenopackets._hpo_ids(ppkt))


def generate(
    omims: list[str],
    out_dir: Path,
    min_hpo: int,
    max_rounds: int,
    *,
    max_hpo: int = 0,
    noise_terms: int = 0,
    seed: int = 0,
) -> dict[str, Path]:
    """create + add-genes per disease, resampling any profile outside the size window,
    then inject noise. Returns {omim: phenopacket path} written into out_dir."""
    out_dir.mkdir(parents=True, exist_ok=True)
    hi = max_hpo or 10**9

    def in_window(n: int) -> bool:
        return min_hpo <= n <= hi

    best: dict[str, tuple[int, Path]] = {}  # omim -> (n_hpo, path); best = closest in-window
    remaining = list(omims)

    def score(n: int) -> tuple[int, int]:
        # prefer in-window; else minimize distance to the window
        dist = 0 if in_window(n) else (min_hpo - n if n < min_hpo else n - hi)
        return (0 if in_window(n) else 1, dist)

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
                accept = in_window(n)
                if accept or omim not in best or score(n) < score(best[omim][0]):
                    dst.write_text(src.read_text())
                    best[omim] = (n, dst)
                if not accept:
                    still.append(omim)
            remaining = still

    off = {o: n for o, (n, _) in best.items() if not in_window(n)}
    if off:
        print(f"warning: outside [{min_hpo},{max_hpo or '∞'}] after {max_rounds} rounds "
              f"(kept closest): {off}")

    if noise_terms > 0:
        _add_noise.universe = sorted({t for ts in _hpoa_disease_terms().values() for t in ts})
        _add_noise.labels = _hp_labels()
        _add_noise.disease_terms = _hpoa_disease_terms()
        for omim, (_, p) in best.items():
            k = _add_noise(p, omim, noise_terms, seed)
            if k < noise_terms:
                print(f"warning: only added {k}/{noise_terms} noise terms to {omim}")

    return {o: p for o, (_, p) in best.items()}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--corpus",
        choices=["synthetic", "synthetic_hard"],
        default=None,
        help="output corpus (default: synthetic, or synthetic_hard with --hard)",
    )
    ap.add_argument(
        "--hard",
        action="store_true",
        help="preset: rare/recent + sparse + noisy → synthetic_hard "
        "(min-hpo 3, max-hpo 12, noise-terms 2) so Exomiser's baseline fails",
    )
    ap.add_argument("--n", type=int, default=20, help="number of synthetic cases")
    ap.add_argument("--seed", type=int, default=0, help="disease-selection + noise seed")
    ap.add_argument("--min-hpo", type=int, default=None, help="min terms in the p2p profile (pre-noise)")
    ap.add_argument("--max-hpo", type=int, default=None, help="max terms in the p2p profile, 0=no cap (pre-noise)")
    ap.add_argument("--noise-terms", type=int, default=None, help="off-target distractors injected on top")
    ap.add_argument("--rare", action="store_true", help="rare/recent sparsely-annotated diseases")
    ap.add_argument("--hpoa-min", type=int, default=6, help="rare pool: min HPOA annotations")
    ap.add_argument("--hpoa-max", type=int, default=15, help="rare pool: max HPOA annotations")
    ap.add_argument("--recent-omim", type=int, default=619000, help="rare pool: min OMIM number")
    ap.add_argument("--pool-min-hpo", type=int, default=4, help="default pool HPOA floor")
    ap.add_argument("--extra-min-hpo", type=int, default=15, help="default backfill HPOA floor")
    ap.add_argument("--max-rounds", type=int, default=12, help="resampling rounds")
    args = ap.parse_args()

    hard = args.hard or args.corpus == "synthetic_hard"
    corpus_name = args.corpus or ("synthetic_hard" if args.hard else "synthetic")
    rare = args.rare or hard
    min_hpo = args.min_hpo if args.min_hpo is not None else (3 if hard else 4)
    max_hpo = args.max_hpo if args.max_hpo is not None else (12 if hard else 0)
    noise = args.noise_terms if args.noise_terms is not None else (2 if hard else 0)
    cp = config.corpus(corpus_name)

    needed = [config.HPOA, config.GENES_TO_DISEASE]
    if noise > 0:
        needed.append(config.HP_OBO)
    for p in needed:
        if not p.exists():
            raise SystemExit(f"{p} not found — run ./download.sh first.")
    if not config.PPKTS.exists():
        raise SystemExit(
            f"{config.PPKTS} not found — run ./download.sh (Phenopacket Store) first."
        )

    omims = select_diseases(
        args.n, args.seed, args.pool_min_hpo, args.extra_min_hpo,
        rare=rare, hpoa_min=args.hpoa_min, hpoa_max=args.hpoa_max,
        recent_omim=args.recent_omim,
    )
    print(
        f"[{corpus_name}] selected {len(omims)} diseases "
        f"(rare={rare}, true-hpo=[{min_hpo},{max_hpo or '∞'}], noise={noise}):\n  "
        + "\n  ".join(omims)
    )

    cp.ppkts.mkdir(parents=True, exist_ok=True)
    for stale in cp.ppkts.glob("*.json"):
        stale.unlink()

    written = generate(
        omims, cp.ppkts, min_hpo, args.max_rounds,
        max_hpo=max_hpo, noise_terms=noise, seed=args.seed,
    )
    print(
        f"\nwrote {len(written)} synthetic phenopackets -> {cp.ppkts}/\n"
        f"next: uv run python -m benchmark.prepare --corpus {corpus_name} --min-hpo {min_hpo}"
    )


if __name__ == "__main__":
    main()
