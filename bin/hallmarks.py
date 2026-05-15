#!/usr/bin/env python3

import argparse
import csv
import pysam
from rapidfuzz import fuzz

# === Constantes===
POLYA_MIN_LEN  = 4      # longitud mínima de cola polyA
POLYA_MAX_GAP  = 2      # hueco máximo entre tramos para fusionarlos
POLYA_MIN_GAP  = 3      # separación mínima entre candidatos tras fusión
POLYA_TAIL_PCT = 0.10   # la cola debe estar en el último/primer 10% de la secuencia
L1_MOTIF       = "TTTTTAA"
L1_MOTIF_RC    = "TTAAAAA"  # complemento reverso
L1_FLANK       = 20     # bp a cada lado de la inserción
L1_MIN_RATIO   = 60     # similitud mínima RapidFuzz (%)


def find_polya(sequence):
#Detecta cola polyA o polyT en la secuencia insertada.
#Devuelve True si:
   #Hay un tramo de A o T de >= POLYA_MIN_LEN bp
   #Ese tramo empieza dentro del último o primer POLYA_TAIL_PCT de la secuencia

    seq = sequence.upper()
    n   = len(seq)
    if n == 0:
        return False

    tail_len = max(1, int(n * POLYA_TAIL_PCT))

    for base in ("A", "T"):
        # Encontrar todos los tramos contiguos de la base
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

        # Fusionar candidatos separados por <= POLYA_MAX_GAP bases
        merged = []
        for start, end in candidates:
            if merged and (start - merged[-1][1]) <= POLYA_MAX_GAP:
                merged[-1] = (merged[-1][0], end)
            else:
                merged.append((start, end))

        # Segunda pasada: asegurar separación mínima de POLYA_MIN_GAP
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

        # Cola al final de la secuencia (strand +)
        best_end   = max(final, key=lambda x: x[1])
        length_end = best_end[1] - best_end[0]
        if length_end >= POLYA_MIN_LEN and best_end[0] >= (n - tail_len):
            return True

        # Cola al inicio de la secuencia (strand -)
        best_start   = min(final, key=lambda x: x[0])
        length_start = best_start[1] - best_start[0]
        if length_start >= POLYA_MIN_LEN and best_start[1] <= tail_len:
            return True

    return False


def find_l1_motif(fasta, chrom, pos):
#Busca el motivo consenso de la endonucleasa de L1 (TTTTTAA o RC TTAAAAA)
#en los +/- L1_FLANK bp alrededor de la posición de inserción.
#Usa búsqueda fuzzy con RapidFuzz al L1_MIN_RATIO% de similitud.
#Devuelve True si se encuentra al menos una coincidencia.
    
    start = max(0, pos - L1_FLANK)
    end   = pos + L1_FLANK

    try:
        region = fasta.fetch(chrom, start, end).upper()
    except (ValueError, KeyError):
        # Cromosoma no encontrado o posición fuera de rango
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
        description="Detección de hallmarks de retrotransposición: cola polyA y motivo L1 endonucleasa"
    )
    parser.add_argument("--input",     required=True, help="TSV con columnas: seqId, seqnames, start, vcf_alt")
    parser.add_argument("--reference", required=True, help="FASTA de referencia indexado (hg38 o T2T)")
    parser.add_argument("--output",    required=True, help="TSV de salida: seqId, has_polya, has_l1_motif")
    args = parser.parse_args()

    print("Abriendo referencia...")
    fasta = pysam.FastaFile(args.reference)

    print("Procesando inserciones...")
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

    print(f"Escribiendo resultados ({len(results)} inserciones)...")
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
