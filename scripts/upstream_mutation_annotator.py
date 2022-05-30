"""the script matches the efficiencies of the Ref/ and Alt. Kozak sequences

input: .txt with \t as delimiter:
1 - chromosome (from vcf)
2 - variant position (1-based from vcf)
3 - ID (from vcf)
4 - ref (from vcf)
5 - alt (from vcf)
6 - quality (from vcf)
7 - filter (from vcf)
8 - info (from vcf)
9 - start (0-based, from custom .bed)
10 - end (0-based, from custom .bed)
11 - chain ( + forward, - reverse, from custom .bed)
12 - calculated Kozak position (0 - -6, 1 - -5, 2 - -4, 3 - -3, 4 - -2, 5 - -1, 6 - +1, 7 - +2, 8 - +3, 9 - +4, 10 - +5)
13 - annotation ("upstream", "no_start", 'synonymous'/'missense'/'nonsense', 'Error in annotation')
14 - Kozak type ('AUG_Kozak', 'not_AUG_Kozak')

output: .txt with \t as delimiter:
1 - chromosome (from vcf)
2 - variant position (1-based from vcf)
3 - ID (from vcf)
4 - ref (from vcf)
5 - alt (from vcf)
6 - quality (from vcf)
7 - filter (from vcf)
8 - info (from vcf)
9 - start (0-based, from custom .bed)
10 - end (0-based, from custom .bed)
11 - chain ( + forward, - reverse, from custom .bed)
12 - calculated Kozak position (0 - -6, 1 - -5, 2 - -4, 3 - -3, 4 - -2, 5 - -1, 6 - +1, 7 - +2, 8 - +3, 9 - +4, 10 - +5)
13 - annotation ("upstream", "no_start", 'synonymous'/'missense'/'nonsense', 'Error in annotation')
14 - Kozak type ('AUG_Kozak', 'not_AUG_Kozak')
15 - Ref Kozak efficiency
16 - Ref Kozak efficiency, lower border of confidence interval
17 - Ref Kozak efficiency, upper border of confidence interval
18 - Alt Kozak efficiency
19 - ALt Kozak efficiency, lower border of confidence interval
20 - Alt Kozak efficiency, upper border of confidence interval
21 - change description (getting lower/getting higher/equal)
22 - relative efficiency Eff[alt]/Eff[ref]"""

import argparse
from dna_rna_functions import dna_to_rna_convert

parser = argparse.ArgumentParser()
parser.add_argument("INPUT_PATH", type=str,
                    help="type the path for input txt-file")
parser.add_argument("OUTPUT_PATH", type=str,
                    help="type the path for output txt-file")
parser.add_argument("KOZAK_SEQ_PATH", type=str,
                    help="type the path for the fasta-file, containing Kozak sequences from genome assembly")
parser.add_argument("KOZAK_EFFICIENCY_PATH", type=str,
                    help="type the path for the txt-file with Kozak efficiency data")
args = parser.parse_args()
input_path = args.INPUT_PATH
output_path = args.OUTPUT_PATH
kozak_sequence_path = args.KOZAK_SEQ_PATH
kozak_efficiency_path = args.KOZAK_EFFICIENCY_PATH

kozak_efficiency_dict = {}
kozak_seq_dict = {}

# download the efficiency data as a dictionary (key - Kozak sequence, value - [efficiency, lower border, upper border])
with open(kozak_efficiency_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip().split(sep='\t')
        if '#' not in input_line[0]:
            kozak_efficiency_dict[input_line[0]] = [input_line[1], input_line[2], input_line[3]]

# download the 11-nt Kozak sequences
with open(kozak_sequence_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip()
        if '>' in input_line:
            header = input_line
        else:
            kozak_seq_dict[header] = input_line

# annotation
with open(input_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip().split(sep='\t')
        res = []
        if input_line[13] == 'not_AUG_Kozak':
            # we can not estimate not-AUG-Kozak's, but we write . to maintain the structure
            res = ['.', '.', '.', '.', '.', '.', '.', '.']
            print('\t'.join(input_line))
            print('There is the sequence without classic AUG start codon, I can not work with it')
        else:
            if input_line[3] not in 'ATGC' or input_line[4] not in 'ATGC':
                # we can not estimate the sequences with not-ATGC letters or other symbols,
                # but we write . to maintain the structure
                res = ['.', '.', '.', '.', '.', '.', '.', '.']
                print('\t'.join(input_line))
                print('There is the sequence with strange not-ATGC letters, I can not work with it')
            else:
                if input_line[12] not in ['no_start', 'Error in annotation']:
                    key_kozak = f'>{input_line[0]}:{input_line[8]}-{input_line[9]}'
                    key_kozak_revcom = f'>{input_line[0]}:{input_line[8]}-{input_line[9]}_revcom'

                    if key_kozak in kozak_seq_dict.keys():
                        ref_kozak = dna_to_rna_convert(kozak_seq_dict[key_kozak])
                        if int(input_line[11]) == 0:
                            alt_kozak = dna_to_rna_convert(input_line[4]) + ref_kozak[1:len(ref_kozak)]
                        elif int(input_line[11]) == 10:
                            alt_kozak = ref_kozak[0:len(ref_kozak) - 1] + dna_to_rna_convert(input_line[4])
                        elif int(input_line[11]) in range(6, 9):
                            alt_kozak = None
                        else:
                            alt_kozak = ref_kozak[0:int(input_line[11])] + dna_to_rna_convert(input_line[4]) + \
                                        ref_kozak[int(input_line[11]) + 1:len(ref_kozak)]

                        if ref_kozak is not None and alt_kozak is not None:
                            res = []
                            res.extend(kozak_efficiency_dict[ref_kozak])
                            res.extend(kozak_efficiency_dict[alt_kozak])
                            if int(kozak_efficiency_dict[ref_kozak][0]) > int(kozak_efficiency_dict[alt_kozak][0]):
                                res.append('getting lower')
                            elif int(kozak_efficiency_dict[ref_kozak][0]) < int(kozak_efficiency_dict[alt_kozak][0]):
                                res.append('getting higher')
                            else:
                                res.append('equal')
                            relative_eff = int(kozak_efficiency_dict[alt_kozak][0]) / int(
                                kozak_efficiency_dict[ref_kozak][0])
                            res.append(str(relative_eff))
                        else:
                            # we write . to maintain the structure
                            res = ['.', '.', '.', '.', '.', '.', '.', '.']
                            print('\t'.join(input_line))
                            print('Something went wrong with reconstitution of reference and alternative sequences')
                    elif key_kozak_revcom in kozak_seq_dict.keys():
                        ref_kozak = dna_to_rna_convert(kozak_seq_dict[key_kozak_revcom])
                        if int(input_line[11]) == 0:
                            alt_kozak = dna_to_rna_convert(input_line[4]) + ref_kozak[1:len(ref_kozak)]
                        elif int(input_line[11]) == 10:
                            alt_kozak = ref_kozak[0:len(ref_kozak) - 1] + dna_to_rna_convert(input_line[4])
                        elif int(input_line[11]) in range(6, 9):
                            alt_kozak = None
                        else:
                            alt_kozak = ref_kozak[0:int(input_line[11])] + dna_to_rna_convert(input_line[4]) + \
                                        ref_kozak[int(input_line[11]) + 1:len(ref_kozak)]

                        if ref_kozak is not None and alt_kozak is not None:
                            res = []
                            res.extend(kozak_efficiency_dict[ref_kozak])
                            res.extend(kozak_efficiency_dict[alt_kozak])
                            if kozak_efficiency_dict[ref_kozak][0] > kozak_efficiency_dict[alt_kozak][0]:
                                res.append('getting lower')
                            elif kozak_efficiency_dict[ref_kozak][0] < kozak_efficiency_dict[alt_kozak][0]:
                                res.append('getting higher')
                            else:
                                res.append('equal')
                            relative_eff = int(kozak_efficiency_dict[alt_kozak][0]) / int(
                                kozak_efficiency_dict[ref_kozak][0])
                            res.append(str(relative_eff))
                        else:
                            # we write . to maintain the structure
                            res = ['.', '.', '.', '.', '.', '.', '.', '.']
                            print('\t'.join(input_line))
                            print('Something went wrong with reconstitution of reference and alternative sequences')

                    else:
                        # we write . to maintain the structure
                        res = ['.', '.', '.', '.', '.', '.', '.', '.']
                else:
                    # we write . to maintain the structure
                    res = ['.', '.', '.', '.', '.', '.', '.', '.']

        input_line.extend(res)
        with open(output_path, 'a') as output_file:
            output_file.write('\t'.join(input_line) + '\n')
