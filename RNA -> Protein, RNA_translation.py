def translate_rna(rna):

    codons = {

        "UUU": ("F", "Phe", "Fenylalanin"),
        "UUC": ("F", "Phe", "Fenylalanin"),

        "UUA": ("L", "Leu", "Leucin"),
        "UUG": ("L", "Leu", "Leucin"),
        "CUU": ("L", "Leu", "Leucin"),
        "CUC": ("L", "Leu", "Leucin"),
        "CUA": ("L", "Leu", "Leucin"),
        "CUG": ("L", "Leu", "Leucin"),

        "AUU": ("I", "Ile", "Isoleucin"),
        "AUC": ("I", "Ile", "Isoleucin"),
        "AUA": ("I", "Ile", "Isoleucin"),

        "AUG": ("M", "Met", "Metionin"),

        "GUU": ("V", "Val", "Valin"),
        "GUC": ("V", "Val", "Valin"),
        "GUA": ("V", "Val", "Valin"),
        "GUG": ("V", "Val", "Valin"),

        "UCU": ("S", "Ser", "Serin"),
        "UCC": ("S", "Ser", "Serin"),
        "UCA": ("S", "Ser", "Serin"),
        "UCG": ("S", "Ser", "Serin"),
        "AGU": ("S", "Ser", "Serin"),
        "AGC": ("S", "Ser", "Serin"),

        "CCU": ("P", "Pro", "Prolin"),
        "CCC": ("P", "Pro", "Prolin"),
        "CCA": ("P", "Pro", "Prolin"),
        "CCG": ("P", "Pro", "Prolin"),

        "ACU": ("T", "Thr", "Treonin"),
        "ACC": ("T", "Thr", "Treonin"),
        "ACA": ("T", "Thr", "Treonin"),
        "ACG": ("T", "Thr", "Treonin"),

        "GCU": ("A", "Ala", "Alanin"),
        "GCC": ("A", "Ala", "Alanin"),
        "GCA": ("A", "Ala", "Alanin"),
        "GCG": ("A", "Ala", "Alanin"),

        "UAU": ("Y", "Tyr", "Tyrosin"),
        "UAC": ("Y", "Tyr", "Tyrosin"),

        "CAU": ("H", "His", "Histidin"),
        "CAC": ("H", "His", "Histidin"),

        "CAA": ("Q", "Gln", "Glutamin"),
        "CAG": ("Q", "Gln", "Glutamin"),

        "AAU": ("N", "Asn", "Asparagin"),
        "AAC": ("N", "Asn", "Asparagin"),

        "AAA": ("K", "Lys", "Lysin"),
        "AAG": ("K", "Lys", "Lysin"),

        "GAU": ("D", "Asp", "Asparaginsyra"),
        "GAC": ("D", "Asp", "Asparaginsyra"),

        "GAA": ("E", "Glu", "Glutaminsyra"),
        "GAG": ("E", "Glu", "Glutaminsyra"),

        "UGU": ("C", "Cys", "Cystein"),
        "UGC": ("C", "Cys", "Cystein"),

        "UGG": ("W", "Trp", "Tryptofan"),

        "CGU": ("R", "Arg", "Arginin"),
        "CGC": ("R", "Arg", "Arginin"),
        "CGA": ("R", "Arg", "Arginin"),
        "CGG": ("R", "Arg", "Arginin"),
        "AGA": ("R", "Arg", "Arginin"),
        "AGG": ("R", "Arg", "Arginin"),

        "GGU": ("G", "Gly", "Glycin"),
        "GGC": ("G", "Gly", "Glycin"),
        "GGA": ("G", "Gly", "Glycin"),
        "GGG": ("G", "Gly", "Glycin"),

        "UAA": ("STOP", "STOP", "STOP"),
        "UAG": ("STOP", "STOP", "STOP"),
        "UGA": ("STOP", "STOP", "STOP")
    }

    protein_1 = []
    protein_3 = []
    protein_full = []

    print("RNA-sekvens:", rna)
    print()

    for i in range(0, len(rna), 3):

        codon = rna[i:i+3]

        if len(codon) < 3:
            break

        if codon in codons:

            one, three, full = codons[codon]

            print(f"{codon} → {one} ({three}) - {full}")

            if one == "STOP":
                break

            protein_1.append(one)
            protein_3.append(three)
            protein_full.append(full)

        else:
            print(codon, "→ Okänt kodon")

    print()
    print("Protein (1 bokstav):", "".join(protein_1))
    print("Protein (3 bokstäver):", "-".join(protein_3))
    print("Protein (hela namn):", " - ".join(protein_full))


rna = "AUGGCUUUUGGAGAAUAA"

translate_rna(rna)
