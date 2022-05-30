""" add to the file 2 columns with CLNSIG and gene ID (from the 8th column with INFO)

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
13 - Kozak type ('AUG_Kozak', 'not_AUG_Kozak')
14 - annotation ("upstream", "no_start", 'synonymous'/'missense'/'nonsense', 'Error in annotation')
15 - Ref Kozak efficiency
16 - Ref Kozak efficiency, lower border of confidence interval
17 - Ref Kozak efficiency, upper border of confidence interval
18 - Alt Kozak efficiency
19 - ALt Kozak efficiency, lower border of confidence interval
20 - Alt Kozak efficiency, upper border of confidence interval
21 - change description (getting lower/getting higher/equal)
22 - relative efficiency Eff[alt]/Eff[ref]

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
13 - Kozak type ('AUG_Kozak', 'not_AUG_Kozak')
14 - annotation ("upstream", "no_start", 'synonymous'/'missense'/'nonsense', 'Error in annotation')
15 - Ref Kozak efficiency
16 - Ref Kozak efficiency, lower border of confidence interval
17 - Ref Kozak efficiency, upper border of confidence interval
18 - Alt Kozak efficiency
19 - ALt Kozak efficiency, lower border of confidence interval
20 - Alt Kozak efficiency, upper border of confidence interval
21 - change description (getting lower/getting higher/equal)
22 - relative efficiency Eff[alt]/Eff[ref]
23 - clinical significance
24 - gene ID """

import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("INPUT_PATH", type=str,
                    help="type the path for input txt-file")
parser.add_argument("OUTPUT_PATH", type=str,
                    help="type the path for output txt-file")
args = parser.parse_args()
input_path = args.INPUT_PATH
output_path = args.OUTPUT_PATH

clnsig_pattern = r'CLNSIG=(\w+\/*\w*)'
geneID_pattern = r'GENEINFO=(\w+)'

found_clnsig = ''
found_geneID = ''

with open(input_path, 'r') as input_file:
    for line in input_file.readlines():
        input_line = line.strip().split(sep='\t')
        found_clnsig = ";".join(re.findall(string=input_line[7], pattern=clnsig_pattern))
        found_geneID = ";".join(re.findall(string=input_line[7], pattern=geneID_pattern))
        input_line.extend([found_clnsig, found_geneID])
        with open(output_path, 'a') as output_file:
            output_file.write('\t'.join(input_line) + '\n')