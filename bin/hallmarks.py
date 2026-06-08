#!/usr/bin/env python3

import argparse
import csv
import pysam
from rapidfuzz import fuzz

# Constants
POLYA_MIN_LEN  = 4      # min polyA tail length
POLYA_MAX_GAP  = 2      # max gap between runs to merge
POLYA_MIN_GAP  = 3      # min separation between candidates after merge
POLYA_TAIL_PCT = 0.10   # tail must be within last/first 10% of the sequence
L1_MOTIF       = "TTTTTAA"
L1_MOTIF_RC    = "TTAAAAA"  # reverse complement
L1_FLANK       = 20     # bp flanking the insertion site
L1_MIN_RATIO   = 60     # min RapidFuzz similarity (%)


def find_polya(sequence):
    # Detect polyA or polyT tail in the inserted sequence.
    # Returns True if a run of A or T >= POLYA_MIN_LEN bp
    # starts within the last or first POLYA_TAIL_PCT of the sequence.

    seq = sequence.upper()
    n   = len(seq)
    if n == 0:
        return False

    tail_len = max(1, int(n * POLYA_TAIL_PCT))

    for base in ("A", "T"):
        # Find all contiguous runs of the base
        candidates = []
        i = 0
        while i < n:
            if seq[i] == base:
                j = i
                while j < n and seq[j] == base:
                    j += 1
                candidates.append((i, j))
                i = j
            else:
                i += 1

        # Merge runs separated by <= POLYA_MAX_GAP bases
        merged = []
        for start, end in candidates:
            if merged and (start - merged[-1][1]) <= POLYA_MAX_GAP:
                merged[-1] = (merged[-1][0], end)
            else:
                merged.append((start, end))

        # Second pass: enforce minimum gap between candidates
        final = []
        for idx, (start, end) in enumerate(merged):
            if idx == 0:
                final.append((start, end))
                continue
            prev_end = final[-1][1]
            if (start - prev_end) >= POLYA_MIN_GAP:
                final.append((start, end))
            else:
                final[-1] = (final[-1][0], end)

        if not final:
            continue

        # Tail at the end of the sequence (+ strand)
        best_end   = max(final, key=lambda x: x[1])
        length_end = best_end[1] - best_end[0]
        if length_end >= POLYA_MIN_LEN and best_end[0] >= (n - tail_len):
            return True

        # Tail at the start of the sequence (- strand)
        best_start   = min(final, key=lambda x: x[0])
        length_start = best_start[1] - best_start[0]
        if length_start >= POLYA_MIN_LEN and best_start[1] <= tail_len:
            return True

    return False


def find_l1_motif(fasta, chrom, pos):
    # Search for the L1 endonuclease consensus motif (TTTTTAA / RC TTAAAAA)
    # in +/- L1_FLANK bp around the insertion site.
    # Uses fuzzy matching via RapidFuzz at L1_MIN_RATIO% similarity.
    # Returns True if at least one match is found.

    start = max(0, pos - L1_FLANK)
    end   = pos + L1_FLANK

    try:
        region = fasta.fetch(chrom, start, end).upper()
    except (ValueError, KeyError):
        return False

    if not region:
        return False

    motif_len = len(L1_MOTIF)

    for i in range(len(region) - motif_len + 1):
        window = region[i:i + motif_len]
        if fuzz.ratio(window, L1_MOTIF) >= L1_MIN_RATIO:
            return True
        if fuzz.ratio(window, L1_MOTIF_RC) >= L1_MIN_RATIO:
            return True

    return False


def main():
    parser = argparse.ArgumentParser(
        description="Detect TPRT retrotransposition hallmarks: polyA tail and L1 endonuclease motif"
    )
    parser.add_argument("--input",     required=True, help="TSV with columns: seqId, seqnames, start, vcf_alt")
    parser.add_argument("--reference", required=True, help="Indexed reference FASTA (hg38 or T2T)")
    parser.add_argument("--output",    required=True, help="Output TSV: seqId, has_polya, has_l1_motif")
    args = parser.parse_args()

    print("Loading reference...", flush=True)
    fasta = pysam.FastaFile(args.reference)

    print("Processing insertions...", flush=True)
    results = []

    with open(args.input, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            seq_id = row["seqId"]
            chrom  = row["seqnames"]
            pos    = int(row["start"])
            seq    = row["vcf_alt"]

            has_polya    = find_polya(seq)
            has_l1_motif = find_l1_motif(fasta, chrom, pos)

            results.append({
                "seqId":        seq_id,
                "has_polya":    str(has_polya),
                "has_l1_motif": str(has_l1_motif)
            })

    fasta.close()

    print(f"Done, {len(results)} insertions processed", flush=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["seqId", "has_polya", "has_l1_motif"],
            delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(results)


if __name__ == "__main__":
    main()
