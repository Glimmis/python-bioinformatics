def base_counter(dna):
    print("DNA-sekvens:", dna)
    print("Totalt antal baser =", len(dna))
    print()

    for bas in "ATCG":
        print(bas, "=", dna.count(bas))

dna = "ATCGGCTAGCGCGATAT"

base_counter(dna)
