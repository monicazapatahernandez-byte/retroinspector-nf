#! /usr/bin/env python3
import argparse
import subprocess
import tempfile
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_bed")
    parser.add_argument("sva_types")
    parser.add_argument("output_bed")
    args = parser.parse_args()

    tmpFasta = args.output_bed + ".tmp.fasta"
    svafDb = args.output_bed + ".blastdb"
    blastOut = args.output_bed + ".blast"

    # Read rm.bed and save SVA-F in a temp fasta file
    with open(args.input_bed) as inputBed, \
            open(tmpFasta, "w") as fasta:
        fastaLines = []
        for line in inputBed.readlines():
            line = line.rstrip().split("\t")
            name = line[3]
            seqId = line[6]
            sequence = line[17]
            fastaLines.append(f">{seqId}\n{sequence}\n")
        fasta.writelines(fastaLines)

    # Create blastn db from sva_types and run blastn
    print(f"svaf db: {svafDb}")
    subprocess.run(["makeblastdb", "-in", args.sva_types,
                   "-out", svafDb, "-dbtype", "nucl"])
    subprocess.run(["blastn", "-task", "blastn", "-evalue", "1e-20",
                   "-db", svafDb, "-query", tmpFasta,
                   "-outfmt", "6", "-out", blastOut])

    # Create dict from blastn results
    blastResults = {}
    with open(blastOut) as blast:
        for line in blast.readlines():
            line = line.split("\t")
            seqId = line[0]
            if seqId not in blastResults:
                blastResults[seqId] = line[1]

    # Read rm.bed again, replace SVA-F with SVA-F1 if in blast results
    with open(args.input_bed) as inputBed, \
        open(args.output_bed, "w") as o:
        lines = []
        for line in inputBed.readlines():
            line = line.rstrip().split("\t")
            seqId = line[6]
            if seqId in blastResults:
                line[3] = blastResults[seqId]
            lines.append("\t".join(line) + "\n")
        o.writelines(lines)

if __name__ == "__main__":
    main()
