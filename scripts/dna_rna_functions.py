# to convert DNA sequence in reverse complement DNA
def dna_rev_com(seq):
    seq_to_revcom = reversed(seq.upper())
    res_seq = ""
    com_dict = {"A": "T", "T": "A", "G": "C", "C": "G"}
    for letter in seq_to_revcom:
        res_seq += com_dict[letter]
    return res_seq


# to convert DNA sequence in RNA
def dna_to_rna_convert(seq):
    if len(set(seq)) > 4:
        raise ValueError("Your sequence seem to have more than A, T, G, C")
    res_seq = ""
    for letter in seq.upper():
        if letter == "T":
            res_seq += "U"
        else:
            res_seq += letter
    return res_seq


# to annotate mutation in a protein-coding sequence
def mutation_type(ref_triplet, alt_triplet):
    genetic_code = {'UUU': "F", 'UUC': "F", 'UUA': "L", 'UUG': "L",
                    'UCU': "S", 'UCC': "S", 'UCA': "S", 'UCG': "S",
                    'UAU': "Y", 'UAC': 'Y', 'UAA': 'stop', 'UAG': 'stop',
                    'UGU': 'C', 'UGC': 'C', 'UGA': 'stop', 'UGG': 'W',
                    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
                    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
                    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
                    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
                    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
                    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
                    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
                    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
    if genetic_code[alt_triplet] == genetic_code[ref_triplet]:
        current_mutation_type = 'synonymous'
    elif genetic_code[alt_triplet] == 'stop':
        current_mutation_type = 'nonsense'
    else:
        current_mutation_type = 'missense'
    return current_mutation_type
