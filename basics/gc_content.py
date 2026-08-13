#!/usr/bin/env python3
"""Calculate the GC content of a DNA sequence, as a percentage."""

import argparse


def gc_content(dna):
    if not dna:
        return 0.0
    return (dna.count("G") + dna.count("C")) / len(dna) * 100


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("dna", nargs="?", default="ATCGGCTAGCGCGATAT",
                        help="DNA sequence (default: a short example)")
    args = parser.parse_args()
    dna = args.dna.upper()

    print("DNA sequence:", dna)
    print("GC content: {:.2f} %".format(gc_content(dna)))


if __name__ == "__main__":
    main()
