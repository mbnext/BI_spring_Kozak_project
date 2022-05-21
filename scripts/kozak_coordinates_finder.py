# input: genome annotation .gff (without mitochondrial DNA, containns only transcript and corresponding exons and CDS
# script finds the first exon of evey transcript, then find start of its CDS and make the interval of Kozak sequence [-6; +5]
# if it is + chain, the interval is [-7, +5) (0 based, right border is not included)
# if it is - chain, the interval is [-6, +6) (0 based, right border is not included)
# output: .bed file: $1 - chromosome, $2 - start, $3 - end, $4 - chain orientation

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("INPUT_PATH", type=str,
                    help="type the path for input file")
parser.add_argument("OUTPUT_PATH", type=str,
                    help="type the path for output file")
args = parser.parse_args()
input_path = args.INPUT_PATH
output_path = args.OUTPUT_PATH

is_transcript = False
is_exon = False
output_res = []
with open (input_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip().split(sep='\t')
        if is_transcript and is_exon:
            if 'CDS' in input_line[2]:
                is_transcript = False
                is_exon = False
                name = input_line[0]
                if "+" in input_line[6]:
                    res = [name, int(input_line[3]) - 7, int(input_line[3]) + 5, input_line[6]]
                    # plus-chain (right border is not included)
                else:
                    res = [name, int(input_line[4]) - 6, int(input_line[4]) + 6, input_line[6]]
                    # minus-chain (right border is not included)
                if res not in output_res:
                    # save only unique
                    output_res.append(res)
        if not is_transcript and 'transcript' in input_line[2]:
            is_transcript = True
        if is_transcript and 'exon' in input_line[2]:
            is_exon = True

with open (output_path, 'a') as output_file:
    for line in output_res:
        output_file.write('\t'.join(str(el) for el in line) + '\n')