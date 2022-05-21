# script for "crude" annotation:
# the position is in range from 0 to 5 = "upstream"
# the position is in range from 6 to 8, in the start codon = "no_start"
# the position is 9 or 10 and Kozak sequences plus 1 nt are used to understand the change in 2nd codon = 'synonymous'/'missense'/'nonsense'
# if the Ref letter from the .vcf file is not on the dedicated position in the Kozak sequence extracted from the genome = 'Error in annotation'
# the last info is written in the file with errors too.

# input: .txt with \t as delimiter:
# 1 - chromosome (from vcf)
# 2 - variant position (1-based from vcf)
# 3 - ID (from vcf)
# 4 - ref (from vcf)
# 5 - alt (from vcf)
# 6 - quality (from vcf)
# 7 - filter (from vcf)
# 8 - info (from vcf)
# 9 - start (0-based, from custom .bed)
# 10 - end (0-based, from custom .bed)
# 11 - chain ( + forward, - reverse, from custom .bed)
# 12 - calculated Kozak position (0 - -6, 1 - -5, 2 - -4, 3 - -3, 4 - -2, 5 - -1, 6 - +1, 7 - +2, 8 - +3, 9 - +4, 10 - +5)

# output: .txt with \t as delimiter:
# 1 - chromosome (from vcf)
# 2 - variant position (1-based from vcf)
# 3 - ID (from vcf)
# 4 - ref (from vcf)
# 5 - alt (from vcf)
# 6 - quality (from vcf)
# 7 - filter (from vcf)
# 8 - info (from vcf)
# 9 - start (0-based, from custom .bed)
# 10 - end (0-based, from custom .bed)
# 11 - chain ( + forward, - reverse, from custom .bed)
# 12 - calculated Kozak position (0 - -6, 1 - -5, 2 - -4, 3 - -3, 4 - -2, 5 - -1, 6 - +1, 7 - +2, 8 - +3, 9 - +4, 10 - +5)
# 13 - annotation ("upstream", "no_start", 'synonymous'/'missense'/'nonsense', 'Error in annotation')
# 14 - Kozak type ('AUG_Kozak', 'not_AUG_Kozak')

import argparse
from dna_rna_functions import mutation_type
from dna_rna_functions import dna_to_rna_convert
from dna_rna_functions import dna_rev_com

parser = argparse.ArgumentParser()
parser.add_argument("INPUT_PATH", type=str,
                    help="type the path for input txt-file")
parser.add_argument("OUTPUT_PATH", type=str,
                    help="type the path for output txt-file")
parser.add_argument("KOZAK_SEQ_PATH", type=str,
                    help="type the path for the fasta-file, containing Kozak sequences +1 nt from genome assembly")
parser.add_argument("KOZAK_COORDINATES_PATH", type=str,
                    help="type the path for the bed-file, containing Kozak bed coordinates")
parser.add_argument("ERROR_PATH", type=str,
                    help="type the path for the txt-file for error log")
args = parser.parse_args()
input_path = args.INPUT_PATH
output_path = args.OUTPUT_PATH
kozak_plus_1_fasta_path = args.KOZAK_SEQ_PATH
kozak_coordinates_path = args.KOZAK_COORDINATES_PATH
error_path = args.ERROR_PATH

kozak_plus_1_dict = {}
kozak_without_ATG = {}
kozak_coordinates = {}
kozak_plus_1_modified = {}

# download the 12-nt "Kozak" sequences
with open(kozak_plus_1_fasta_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip()
        if '>' in input_line:
            header = input_line
        else:
            if 'N' not in input_line:
                kozak_plus_1_dict[header] = input_line

# download the coordinates from .bed file in a dictionary with the keys corresponding to .fasta names
with open(kozak_coordinates_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip().split(sep='\t')
        key_coord = f'>{input_line[0]}:{input_line[1]}-{input_line[2]}'
        kozak_coordinates[key_coord] = input_line[3]

# classify the Kozak's (chain + or -, contains AUG or do not contain AUG
for key in kozak_plus_1_dict.keys():
    val = kozak_plus_1_dict[key]
    if key in kozak_coordinates.keys():
        if kozak_coordinates[key] == '+':
            if val[6:9] == 'ATG':
                kozak_plus_1_modified[key] = val
            else:
                kozak_plus_1_modified[f'{key}_noAUG'] = val
                kozak_without_ATG[key] = val
        else:
            if val[3:6] == 'CAT':
                kozak_plus_1_modified[f'{key}_revcom'] = dna_rev_com(val)
            else:
                kozak_plus_1_modified[f'{key}_revcom_noAUG'] = dna_rev_com(val)
                kozak_without_ATG[f'{key}_revcom'] = dna_rev_com(val)

# write the file with not_AUG_Kozak
if len(kozak_without_ATG.keys()) > 0:
    print("I have found some Kozak sequences without standard AUG start codon, I am writing them to the file "
          "but they are used in the analysis")
    with open("found_Kozak_without_ATG.fasta", 'w') as error1_file:
        for key in kozak_without_ATG.keys():
            error1_file.write(key + '\n')
            error1_file.write(kozak_without_ATG[key] + '\n')


# annotation function
def crude_annotator(input_string_as_list, kozak_plus_one_dictionary, key):
    inp_data = input_string_as_list
    res1 = ''
    res2 = ''
    ref = ''
    alt = ''
    position = int(inp_data[11])
    ref_letter_at_position = input_line[3]
    alt_letter_at_position = input_line[4]

    if key in kozak_plus_one_dictionary.keys():
        if '_noAUG' in key:
            res2 = 'not_AUG_Kozak'
        else:
            res2 = 'AUG_Kozak'
    else:
        res2 = 'This Kozak sequence not found in the provided sequences'

    if int(inp_data[11]) < 6:
        res1 = 'upstream'
    elif int(inp_data[11]) in range(6, 9):
        res1 = 'no_start'
    else:
        if key in kozak_plus_one_dictionary.keys():
            ref = kozak_plus_one_dictionary[key][9:12]
            if ref_letter_at_position in 'ATGC' and alt_letter_at_position in 'ATGC':
                if position == 9:
                    if ref_letter_at_position != ref[0]:
                        print(input_line[8], input_line[9], ref_letter_at_position, ref[0])
                        print("Reference base from .vcf does not match to the extracted reference sequence")
                        with open(error_path, 'a') as error_file:
                            error_file.write(
                                '\t'.join([input_line[8], input_line[9], input_line[3], ref[0]]) + '\n')
                        alt = None
                    else:
                        alt = alt_letter_at_position + ref[1:len(ref)]
                elif position == 10:
                    if ref_letter_at_position != ref[1]:
                        print(input_line[8], input_line[9], ref_letter_at_position, ref[1])
                        print("Reference base from .vcf does not match to the extracted reference sequence")
                        with open(error_path, 'a') as error_file:
                            error_file.write(
                                '\t'.join([input_line[8], input_line[9], input_line[3], ref[0]]) + '\n')
                        alt = None
                    else:
                        alt = ref[0] + alt_letter_at_position + ref[-1]
                    # последняя буква триплета в последовательность Козак не входит, ее не проверяем
            else:
                ref, alt = None, None

            if ref is not None and alt is not None:
                res1 = mutation_type(alt_triplet=dna_to_rna_convert(alt), ref_triplet=dna_to_rna_convert(ref))
            else:
                res1 = 'Error in annotation'

    return res1, res2


# the final analysis
with open(input_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip().split(sep='\t')
        key_kozak = f'>{input_line[0]}:{input_line[8]}-{input_line[9]}'
        key_kozak_revcom = f'>{input_line[0]}:{input_line[8]}-{input_line[9]}_revcom'
        key_kozak_noAUG = f'>{input_line[0]}:{input_line[8]}-{input_line[9]}_noAUG'
        key_kozak_noAUG_revcom = f'>{input_line[0]}:{input_line[8]}-{input_line[9]}_revcom_noAUG'

        if key_kozak in kozak_plus_1_modified.keys():
            res1, res2 = crude_annotator(input_line, kozak_plus_1_modified, key_kozak)
        elif key_kozak_revcom in kozak_plus_1_modified.keys():
            res1, res2 = crude_annotator(input_line, kozak_plus_1_modified, key_kozak_revcom)
        elif key_kozak_noAUG in kozak_plus_1_modified.keys():
            res1, res2 = crude_annotator(input_line, kozak_plus_1_modified, key_kozak_noAUG)
        elif key_kozak_noAUG_revcom in kozak_plus_1_modified.keys():
            res1, res2 = crude_annotator(input_line, kozak_plus_1_modified, key_kozak_noAUG_revcom)
        else:
            res1, res2 = 'Error in variant annotation', 'Error in AUG/not_AUG estimation'

        input_line.append(res1)
        input_line.append(res2)

        with open(output_path, 'a') as output_file:
            output_file.write('\t'.join(input_line) + '\n')