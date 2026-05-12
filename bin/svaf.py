#! /usr/bin/env python
import subprocess
import sys
import os

input_bed = sys.argv[1]
svaf_types = sys.argv[2]
output_bed = sys.argv[3]

output_fasta = "svaf_tmp.fasta"
output_blast = "svaf_tmp.blast6"
db_prefix = "svaf_db"

with open(input_bed) as inputBed, open(output_fasta, "w") as fasta:
    fastaLines = []
    for line in inputBed.readlines():
        line = line.rstrip().split("\t")
        name = line[3]
        seqId = line[6]
        sequence = line[17]
        fastaLines.append(f">{seqId}\n{sequence}\n")
    fasta.writelines(fastaLines)

subprocess.run(["makeblastdb", "-in", svaf_types, "-out", db_prefix, "-dbtype", "nucl"])
subprocess.run(["blastn", "-task", "blastn", "-evalue", "1e-20", "-db", db_prefix,
               "-query", output_fasta, "-outfmt", "6", "-out", output_blast])

blastResults = {}
if os.path.exists(output_blast):
    with open(output_blast) as blast:
        for line in blast.readlines():
            line = line.split("\t")
            seqId = line[0]
            if seqId not in blastResults:
                blastResults[seqId] = line[1]

with open(input_bed) as inputBed, open(output_bed, "w") as o:
    lines = []
    for line in inputBed.readlines():
        line = line.rstrip().split("\t")
        seqId = line[6]
        if seqId in blastResults:
            line[3] = blastResults[seqId]
        lines.append("\t".join(line) + "\n")
    o.writelines(lines)
