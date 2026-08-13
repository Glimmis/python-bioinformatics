#!/usr/bin/env python3
"""
Generate the synthetic example contig used in the README.

The contig is SYNTHETIC, not a real genome: random intergenic sequence with a
handful of protein-coding ORFs planted at known positions on both strands.
Because the sequence is generated from a fixed random seed, running this script
reproduces example_contig.fasta byte for byte.

Planted genes (forward strand unless noted):
    gene_A  180 aa
    gene_B  240 aa
    gene_C  120 aa   (reverse strand)
    gene_D  310 aa

Usage:
    python make_example_contig.py > example_contig.fasta
"""

import random
import sys

SEED = 20260813
STOPS = {"TAA", "TAG", "TGA"}
SENSE_CODONS = [
    a + b + c
    for a in "ACGT" for b in "ACGT" for c in "ACGT"
    if a + b + c not in STOPS
]

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def random_intergenic(rng, length):
    """Random sequence with any in-frame ATG removed, so no stray ORFs appear."""
    bases = [rng.choice("ACGT") for _ in range(length)]
    sequence = "".join(bases)
    return sequence.replace("ATG", "ATC")


def make_gene(rng, protein_length):
    """ATG + protein_length-1 sense codons + stop codon."""
    body = "".join(rng.choice(SENSE_CODONS) for _ in range(protein_length - 1))
    return "ATG" + body + rng.choice(sorted(STOPS))


def build_contig():
    rng = random.Random(SEED)
    parts = [random_intergenic(rng, 250)]

    for protein_length, on_reverse in [(180, False), (240, False),
                                       (120, True), (310, False)]:
        gene = make_gene(rng, protein_length)
        parts.append(reverse_complement(gene) if on_reverse else gene)
        parts.append(random_intergenic(rng, 200))

    return "".join(parts)


def write_fasta(sequence, header, out=sys.stdout, width=70):
    print(">" + header, file=out)
    for start in range(0, len(sequence), width):
        print(sequence[start:start + width], file=out)


if __name__ == "__main__":
    contig = build_contig()
    write_fasta(
        contig,
        "example_contig synthetic test contig, {} bp, 4 planted ORFs "
        "(3 forward, 1 reverse)".format(len(contig)),
    )
