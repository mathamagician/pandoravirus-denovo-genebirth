"""
BLASTp Complete — Remaining Hypothetical Genes
===============================================
BLASTs all 902 hypothetical P. salinus proteins against NCBI nr,
skipping the 300 already completed in blastp_hypothetical.json.

Estimated time: 6-10 hours wall-clock (602 queries × ~60s each + NCBI latency).
Fully resumable — safe to interrupt and restart at any time.

Usage:
    python scripts/blast_hyp_complete.py

Output:
    results/blastp_hypothetical.json  (extended in-place, same format as before)
    results/blast_hyp_complete_log.txt  (progress log)
"""

import json
import time
import sys
import os
from pathlib import Path
from datetime import datetime

from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML

# ── Config ─────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR     = PROJECT_ROOT / "data"
RESULTS_DIR  = PROJECT_ROOT / "results"

RESULTS_PATH = RESULTS_DIR / "blastp_hypothetical.json"
FASTA_PATH   = DATA_DIR / "proteins_hypothetical.fasta"
LOG_PATH     = RESULTS_DIR / "blast_hyp_complete_log.txt"

Entrez.email  = "Eddie.Bradford@gmail.com"
Entrez.api_key = None   # set your NCBI API key here for 10 req/s

BLAST_EVALUE   = 1e-5
BLAST_MAX_HITS = 10
BLAST_DB       = "nr"
BLAST_DELAY    = 1      # seconds between successful queries


def log(msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    with open(LOG_PATH, "a") as f:
        f.write(line + "\n")


def extract_organism(hit_def):
    import re
    matches = re.findall(r'\[([^\]]+)\]', hit_def)
    return matches[-1] if matches else "Unknown"


def run_remaining():
    # Load all hypothetical proteins from FASTA
    all_proteins = list(SeqIO.parse(str(FASTA_PATH), "fasta"))
    log(f"Total hypothetical proteins in FASTA: {len(all_proteins)}")

    # Load existing results
    if RESULTS_PATH.exists():
        with open(RESULTS_PATH) as f:
            existing = json.load(f)
        log(f"Already completed: {len(existing)}")
    else:
        existing = {}
        log("No existing results — starting fresh")

    # Identify remaining
    done_ids  = set(existing.keys())
    remaining = [p for p in all_proteins if p.id not in done_ids]
    total     = len(all_proteins)
    done      = total - len(remaining)

    log(f"Remaining: {len(remaining)} queries  |  Done: {done}/{total}")

    if not remaining:
        log("All proteins already BLASTed. Nothing to do.")
        return existing

    # BLAST remaining proteins
    for i, protein in enumerate(remaining):
        n = done + i + 1
        log(f"[{n}/{total}] {protein.id}  {protein.description[:70]}")
        t0 = time.time()

        try:
            result_handle = NCBIWWW.qblast(
                "blastp", BLAST_DB, str(protein.seq),
                expect=BLAST_EVALUE,
                hitlist_size=BLAST_MAX_HITS,
                format_type="XML",
            )

            blast_records = NCBIXML.parse(result_handle)
            record = next(blast_records)

            hits = []
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    organism = extract_organism(alignment.hit_def)
                    hits.append({
                        "hit_id":         alignment.hit_id,
                        "hit_def":        alignment.hit_def[:300],
                        "hit_accession":  alignment.accession,
                        "organism":       organism,
                        "subject_length": alignment.length,
                        "evalue":         hsp.expect,
                        "bit_score":      hsp.bits,
                        "score":          hsp.score,
                        "identities":     hsp.identities,
                        "positives":      hsp.positives,
                        "gaps":           hsp.gaps,
                        "align_length":   hsp.align_length,
                        "query_start":    hsp.query_start,
                        "query_end":      hsp.query_end,
                        "pct_identity":   round(hsp.identities / hsp.align_length * 100, 1) if hsp.align_length else 0,
                        "query_coverage": round((hsp.query_end - hsp.query_start + 1) / len(protein.seq) * 100, 1),
                    })
                    break  # top HSP per alignment only

            non_self = [h for h in hits if "pandoravirus" not in h["organism"].lower()]

            existing[protein.id] = {
                "protein_id":           protein.id,
                "description":          protein.description,
                "query_length":         len(protein.seq),
                "n_hits_total":         len(hits),
                "n_hits_non_self":      len(non_self),
                "hits":                 hits,
                "non_self_hits":        non_self,
                "top_organism":         hits[0]["organism"] if hits else None,
                "top_non_self_organism": non_self[0]["organism"] if non_self else None,
            }

            elapsed = time.time() - t0
            top_ns  = non_self[0]["organism"] if non_self else "no non-self hits"
            log(f"  -> {len(hits)} hits ({len(non_self)} non-self) in {elapsed:.1f}s | top non-self: {top_ns}")

        except Exception as e:
            log(f"  ERROR: {e}")
            existing[protein.id] = {
                "protein_id":      protein.id,
                "description":     protein.description,
                "query_length":    len(protein.seq),
                "n_hits_total":    0,
                "n_hits_non_self": 0,
                "hits":            [],
                "non_self_hits":   [],
                "error":           str(e),
            }

        # Save after every query — fully resumable
        with open(RESULTS_PATH, "w") as f:
            json.dump(existing, f, indent=2)

        time.sleep(BLAST_DELAY)

    log(f"COMPLETE. Total entries in results file: {len(existing)}")
    return existing


if __name__ == "__main__":
    log("=" * 60)
    log("blast_hyp_complete.py — starting full hypothetical BLAST")
    log("=" * 60)
    run_remaining()
