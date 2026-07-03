# BLAST annotation add-on for the ORF finder.
#
# Takes the ORFs you already found and asks NCBI whether each predicted protein
# matches a known protein. This is "homology / sequence-similarity search",
# done with BLASTp against NCBI's non-redundant (nr) protein database.
#
# Requirements:
#   pip install biopython
#   An internet connection (queries NCBI's public BLAST server).
#
# NOTE: online BLAST is SLOW (roughly 0.5-3 minutes per sequence) and rate
# limited. Only BLAST the ORFs worth checking (use a length filter) and be
# patient. For hundreds of sequences you would install BLAST locally instead.

import csv
import time

from Bio.Blast import NCBIWWW, NCBIXML


def blast_one_protein(protein_sequence, evalue_cutoff=1e-5, hitlist_size=5):
    """
    Run BLASTp for a single protein sequence against NCBI 'nr'.

    Returns a dict describing the best hit, or None if there was no significant
    hit below evalue_cutoff.
    """
    # Submit the search and wait for NCBI to return XML results.
    result_handle = NCBIWWW.qblast(
        program="blastp",
        database="nr",
        sequence=protein_sequence,
        expect=evalue_cutoff,
        hitlist_size=hitlist_size,
    )

    blast_record = NCBIXML.read(result_handle)
    result_handle.close()

    if len(blast_record.alignments) == 0:
        return None

    # The first alignment is the best-scoring hit.
    best_alignment = blast_record.alignments[0]
    best_hsp = best_alignment.hsps[0]

    if best_hsp.expect > evalue_cutoff:
        return None

    percent_identity = 100.0 * best_hsp.identities / best_hsp.align_length

    return {
        "hit_title": best_alignment.title,
        "e_value": best_hsp.expect,
        "percent_identity": round(percent_identity, 1),
        "align_length": best_hsp.align_length,
    }


def annotate_orfs_with_blast(orfs, min_protein_length=50, evalue_cutoff=1e-5,
                             max_orfs=None, pause_seconds=3):
    """
    Annotate ORFs by BLASTing their proteins against NCBI.

    orfs               : the list of ORF dicts from run_orf_finder / analyze_fasta_file
    min_protein_length : skip ORFs shorter than this (aa) - saves time and noise
    evalue_cutoff      : significance threshold; smaller = stricter
    max_orfs           : optionally cap how many ORFs to BLAST (for a quick test)
    pause_seconds      : delay between requests to respect NCBI rate limits

    Adds keys to each ORF: is_probable_gene, blast_hit, blast_evalue,
    blast_identity. Returns the same list.
    """
    # Only BLAST ORFs that are long enough to be worth checking.
    candidates = [o for o in orfs if o["protein_length_aa"] >= min_protein_length]

    # Longest first - the real genes tend to be the longer ORFs.
    candidates.sort(key=lambda o: o["protein_length_aa"], reverse=True)

    if max_orfs is not None:
        candidates = candidates[:max_orfs]

    print("ORFs to BLAST:", len(candidates),
          "(protein length >=", min_protein_length, "aa)")

    for number, orf in enumerate(candidates, start=1):
        protein = orf["protein_1_letter"]
        print()
        print("[{}/{}] BLASTing {} {}-{} ({} aa)...".format(
            number, len(candidates), orf["strand"],
            orf["start"], orf["end"], orf["protein_length_aa"]))

        try:
            hit = blast_one_protein(protein, evalue_cutoff=evalue_cutoff)
        except Exception as error:
            print("   BLAST failed:", error)
            orf["is_probable_gene"] = None
            orf["blast_hit"] = "ERROR: " + str(error)
            orf["blast_evalue"] = ""
            orf["blast_identity"] = ""
            time.sleep(pause_seconds)
            continue

        if hit is None:
            print("   No significant hit -> probably a random ORF.")
            orf["is_probable_gene"] = False
            orf["blast_hit"] = "No significant hit"
            orf["blast_evalue"] = ""
            orf["blast_identity"] = ""
        else:
            print("   Match:", hit["hit_title"][:70])
            print("   E-value:", hit["e_value"],
                  "| identity:", hit["percent_identity"], "%")
            orf["is_probable_gene"] = True
            orf["blast_hit"] = hit["hit_title"]
            orf["blast_evalue"] = hit["e_value"]
            orf["blast_identity"] = hit["percent_identity"]

        # Be polite to NCBI: wait between requests.
        time.sleep(pause_seconds)

    return orfs


def export_annotated_orfs_to_csv(orfs, filename="orf_annotation.csv"):
    """Export only the ORFs that were BLASTed, with their annotation."""
    annotated = [o for o in orfs if "is_probable_gene" in o]

    if len(annotated) == 0:
        print("No annotated ORFs to export.")
        return

    with open(filename, "w", newline="", encoding="utf-8") as csvfile:
        fieldnames = [
            "sequence_id", "strand", "frame", "start", "end",
            "protein_length_aa", "is_probable_gene",
            "blast_hit", "blast_evalue", "blast_identity",
            "protein_1_letter",
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for orf in annotated:
            writer.writerow(orf)

    print("Annotation exported to", filename)


# ------------------------------------------------------------
# Example usage (run AFTER you have the `orfs` list from orf_finder.py)
# ------------------------------------------------------------
# from orf_finder import run_orf_finder
# orfs = run_orf_finder("phiX174.fasta", min_protein_length=50)
#
# # BLAST the 3 longest ORFs as a quick test (each takes a minute or two):
# annotate_orfs_with_blast(orfs, min_protein_length=50, max_orfs=3)
# export_annotated_orfs_to_csv(orfs, "phiX174_annotation.csv")
