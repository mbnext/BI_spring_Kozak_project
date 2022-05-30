"""if bedtools getfasta extracted the Kozak sequences with 1 additional nt, this script write them as 11-nt Kozak's and
also writes the 'not_AUG_Kozak' in the file."""

import argparse
from dna_rna_functions import dna_rev_com

parser = argparse.ArgumentParser()
parser.add_argument("OUTPUT_PATH", type=str,
                    help="type the path for output txt-file")
parser.add_argument("KOZAK_SEQ_PATH", type=str,
                    help="type the path for the fasta-file, containing Kozak sequences +1 nt from genome assembly")
parser.add_argument("KOZAK_COORDINATES_PATH", type=str,
                    help="type the path for the bed-file, containing Kozak bed coordinates")
args = parser.parse_args()

output_path = args.OUTPUT_PATH
kozak_plus_1_fasta_path = args.KOZAK_SEQ_PATH
kozak_coordinates_path = args.KOZAK_COORDINATES_PATH

kozak_plus_1_dict = {}
kozak_without_ATG = {}
kozak_coordinates = {}
kozak_plus_1_modified = {}
kozak_dict = {}


with open (kozak_plus_1_fasta_path, 'r') as input_file:
    # download the 12-nt 'Kozak' sequences
    for line in input_file.readlines():
        input_line = line.strip()
        if '>' in input_line:
            header = input_line
        else:
            if 'N' not in input_line:
                kozak_plus_1_dict[header] = input_line


with open(kozak_coordinates_path, 'r') as input_file:
    # download the coordinates from .bed file in a dictionary with the keys corresponding to .fasta names
    for line in input_file.readlines():
        input_line = line.strip().split(sep='\t')
        key_coord = f'>{input_line[0]}:{input_line[1]}-{input_line[2]}'
        kozak_coordinates[key_coord] = input_line[3]


for key in kozak_plus_1_dict.keys():
    # classify the Kozak's (chain + or -, contains AUG or do not contain AUG)
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


if len(kozak_without_ATG.keys()) > 0:
    # write the file with not_AUG_Kozak
    print("I have found some Kozak sequences without standard AUG start codon, I am writing them to the file "
          "'found_Kozak_without_ATG.fasta'")
    with open("found_Kozak_without_ATG.fasta", 'w') as error1_file:
        for key in kozak_without_ATG.keys():
            error1_file.write(key + '\n')
            error1_file.write(kozak_without_ATG[key] + '\n')


for key in kozak_plus_1_modified.keys():
    # cut the 12-nt sequence to 11-nt
    kozak_dict[key] = kozak_plus_1_modified[key][0:11]


if len(kozak_dict.keys()) > 0:
    # write the resulted sequences in the file
    with open(output_path, 'w') as output_file:
        for key in kozak_dict.keys():
            output_file.write(key + '\n')
            output_file.write(kozak_dict[key] + '\n')