def gc_halt(dna):
    g = dna.count("G")
    c = dna.count("C")
    gc = g + c
    return gc / len(dna) * 100

dna = "ATCGGCTAGCGCGATAT"

print("DNA-sekvens:", dna)
print("GC-halt:", gc_halt(dna), "%")