#!/usr/bin/env python3
"""
BLAST annotation: functional annotation of predicted ORFs.

Takes the ORF table produced by orf_finder.py and asks NCBI whether each
predicted protein matches a known protein (BLASTp against the non-redundant
'nr' database). This separates probable real genes from ORFs that occur by
chance.

Hits titled "hypothetical / uncharacterized / unnamed protein" are skipped in
favour of a properly named hit further down the list, with the best significant
hit kept as a fallback so a real gene is never lost.

Requirements:
    pip install biopython
    Internet access (queries NCBI's public BLAST server).

Online BLAST is slow (roughly 0.5-3 minutes per sequence) and rate limited.
Use --max-orfs to test on a few sequences first. For hundreds of sequences,
install BLAST+ locally instead.

Usage:
    python blast_annotate.py examples/example_contig_orfs.csv --max-orfs 3
    python blast_annotate.py genome_orfs.csv -o genome_annotation.csv

Also importable:
    from orf_finder import run_orf_finder
    from blast_annotate import annotate_orfs_with_blast, export_annotated_orfs_to_csv
    orfs = run_orf_finder("genome.fasta", min_protein_length=100)
    annotate_orfs_with_blast(orfs, max_orfs=3)
    export_annotated_orfs_to_csv(orfs, "genome_annotation.csv")
"""

import argparse
import csv
import sys
import time

# Title fragments that mean "no useful name".
UNINFORMATIVE_TITLE_MARKERS = (
    "unnamed protein",
    "hypothetical protein",
    "uncharacterized protein",
    "predicted protein",
    "putative uncharacterized",
)

ANNOTATION_CSV_FIELDS = [
    "sequence_id", "strand", "frame", "start", "end", "protein_length_aa",
    "is_probable_gene", "blast_hit", "blast_evalue", "blast_identity",
    "protein_1_letter",
]


def is_informative_title(title):
    """Return True if the hit title looks like a real, named protein."""
    lowered = title.lower()
    return not any(marker in lowered for marker in UNINFORMATIVE_TITLE_MARKERS)


def blast_one_protein(protein_sequence, evalue_cutoff=1e-5, hitlist_size=20):
    """
    BLASTp a single protein against NCBI 'nr'.

    Returns the best significant hit that carries a real name, falling back to
    the best significant hit if every one is uninformative. Returns None if
    there is no significant hit at all.
    """
    from Bio.Blast import NCBIWWW, NCBIXML  # imported lazily so --help works without biopython

    result_handle = NCBIWWW.qblast(
        program="blastp", database="nr", sequence=protein_sequence,
        expect=evalue_cutoff, hitlist_size=hitlist_size,
    )
    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    fallback = None

    for alignment in blast_record.alignments:
        best_hsp = alignment.hsps[0]
        if best_hsp.expect > evalue_cutoff:
            continue

        hit = {
            "hit_title": alignment.title,
            "e_value": best_hsp.expect,
            "percent_identity": round(
                100.0 * best_hsp.identities / best_hsp.align_length, 1),
            "align_length": best_hsp.align_length,
            "named": is_informative_title(alignment.title),
        }

        if fallback is None:
            fallback = hit
        if hit["named"]:
            return hit

    return fallback


def annotate_orfs_with_blast(orfs, min_protein_length=50, evalue_cutoff=1e-5,
                             max_orfs=None, pause_seconds=3):
    """
    BLAST the ORF proteins and record the results on each ORF dict.

    Adds: is_probable_gene, blast_hit, blast_evalue, blast_identity.
    Longest ORFs are processed first, since real genes tend to be longer.
    """
    candidates = [o for o in orfs if int(o["protein_length_aa"]) >= min_protein_length]
    candidates.sort(key=lambda o: int(o["protein_length_aa"]), reverse=True)

    if max_orfs is not None:
        candidates = candidates[:max_orfs]

    print("ORFs to BLAST:", len(candidates),
          "(protein length >=", min_protein_length, "aa)")

    for number, orf in enumerate(candidates, start=1):
        print()
        print("[{}/{}] BLASTing {} {}-{} ({} aa)...".format(
            number, len(candidates), orf["strand"], orf["start"], orf["end"],
            orf["protein_length_aa"]))

        try:
            hit = blast_one_protein(orf["protein_1_letter"],
                                    evalue_cutoff=evalue_cutoff)
        except Exception as error:
            print("   BLAST failed:", error, file=sys.stderr)
            orf.update(is_probable_gene=None, blast_hit="ERROR: " + str(error),
                       blast_evalue="", blast_identity="")
            time.sleep(pause_seconds)
            continue

        if hit is None:
            print("   No significant hit -> probably a random ORF.")
            orf.update(is_probable_gene=False, blast_hit="No significant hit",
                       blast_evalue="", blast_identity="")
        else:
            print("   Match:", hit["hit_title"][:70])
            print("   E-value:", hit["e_value"],
                  "| identity:", hit["percent_identity"], "%")
            orf.update(is_probable_gene=True, blast_hit=hit["hit_title"],
                       blast_evalue=hit["e_value"],
                       blast_identity=hit["percent_identity"])

        time.sleep(pause_seconds)  # be polite to NCBI

    summarize_annotation(orfs)
    return orfs


def summarize_annotation(orfs):
    """Print a signal-vs-noise summary of the BLASTed ORFs."""
    blasted = [o for o in orfs if "is_probable_gene" in o]

    if not blasted:
        print("\nNo ORFs were BLASTed, so there is nothing to summarize.")
        return

    genes = sum(1 for o in blasted if o["is_probable_gene"] is True)
    noise = sum(1 for o in blasted if o["is_probable_gene"] is False)
    errors = sum(1 for o in blasted if o["is_probable_gene"] is None)

    print()
    print("=" * 50)
    print("Annotation summary")
    print("=" * 50)
    print("{} of {} BLASTed ORFs matched a known protein (probable genes)."
          .format(genes, len(blasted)))
    print("No significant hit (probable random ORFs):", noise)
    if errors:
        print("BLAST errors (not classified):", errors)
    print("=" * 50)


def export_annotated_orfs_to_csv(orfs, filename="orf_annotation.csv"):
    """Write only the BLASTed ORFs, with their annotation, to CSV."""
    annotated = [o for o in orfs if "is_probable_gene" in o]

    if not annotated:
        print("No annotated ORFs to export.", file=sys.stderr)
        return

    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=ANNOTATION_CSV_FIELDS,
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(annotated)

    print("Annotation written to", filename)


def load_orfs_from_csv(filepath):
    """Read an ORF table produced by orf_finder.py."""
    with open(filepath, "r", encoding="utf-8") as csvfile:
        return list(csv.DictReader(csvfile))


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Annotate predicted ORFs by BLASTing them against NCBI nr.",
        epilog="Example: python blast_annotate.py examples/example_contig_orfs.csv "
               "--max-orfs 3 -o annotation.csv",
    )
    parser.add_argument("orf_csv", help="ORF CSV produced by orf_finder.py")
    parser.add_argument("-o", "--output", default="orf_annotation.csv",
                        help="output CSV path (default: %(default)s)")
    parser.add_argument("-m", "--min-length", type=int, default=50,
                        help="skip ORFs shorter than this, in aa (default: %(default)s)")
    parser.add_argument("-n", "--max-orfs", type=int, default=None,
                        help="only BLAST this many ORFs (longest first)")
    parser.add_argument("-e", "--evalue", type=float, default=1e-5,
                        help="E-value cutoff (default: %(default)s)")
    parser.add_argument("--pause", type=float, default=3.0,
                        help="seconds to wait between NCBI requests (default: %(default)s)")
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    orfs = load_orfs_from_csv(args.orf_csv)
    if not orfs:
        print("No ORFs found in", args.orf_csv, file=sys.stderr)
        return 1

    annotate_orfs_with_blast(
        orfs, min_protein_length=args.min_length, evalue_cutoff=args.evalue,
        max_orfs=args.max_orfs, pause_seconds=args.pause,
    )
    export_annotated_orfs_to_csv(orfs, args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
