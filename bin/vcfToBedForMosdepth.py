#! /usr/bin/env python3

import argparse
from pysam import VariantFile

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--svtype", default="")
    parser.add_argument("vcf")
    parser.add_argument("bed")
    args = parser.parse_args()

    with open(args.bed, "w") as bed, VariantFile(args.vcf) as vcf:
        for var in vcf.fetch():
            if args.svtype and var.info["SVTYPE"] != args.svtype:
                continue
            bed.write(f"{var.chrom}\t{max(0, var.pos - 61)}\t{var.stop + 60}\t{var.id}\n")

if __name__ == "__main__":
    main()
