#!/usr/bin/env python3
import sys
import os

demix_dir = sys.argv[1]
demix_files = sorted(
    [
        f for f in os.listdir(demix_dir)
        if f.startswith("demix_") and "barcode" in f
    ],
    key=lambda x: int(x.split("barcode")[-1])
)

summary_file = os.path.join(demix_dir, "demix_summary.txt")

# Create summary file if not exists
if not os.path.exists(summary_file):
    with open(summary_file, "w") as f:
        f.write("barcode\tlineage\tabundance\n")

def parse_demix_file(filepath, barcode):
    lineages = []
    abundances = []

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()

            if line.startswith("lineages"):
                parts = line.split("\t", 1)
                if len(parts) == 2 and parts[1].strip():
                    lineages = parts[1].split()
                else:
                    print(f"[WARN] No lineage data in {filepath}")

            elif line.startswith("abundances"):
                parts = line.split("\t", 1)
                if len(parts) == 2 and parts[1].strip():
                    abundances = parts[1].split()
                else:
                    print(f"[WARN] No abundance data in {filepath}")

    # sanity check
    if len(lineages) != len(abundances):
        print(f"[WARN] Lineage/abundance mismatch in barcode {barcode}")
        return None

    return list(zip(lineages, abundances))

with open(summary_file, "a") as out:
    for filename in demix_files:
        barcode = filename.split("barcode")[-1]
        demix_file = os.path.join(demix_dir, filename)

        results = parse_demix_file(demix_file, barcode)

        if results is None:
            continue

        for lineage, abundance in results:
            out.write(f"{barcode}\t{lineage}\t{abundance}\n")
