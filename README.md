# Prediction of pathogenicity of genetic variants in Kozak sequences

Student: Marianna Baranovskaia

Supervisors: Yury Barbitov (Bioinformatics Institute), Michail Skoblov (Research Center of Medical Genetics)

## Introduction
Kozak sequence is a consensus nucleotide environment of the start codon in the most of the eukaryotic mRNAs, involved in the translation initiation. Marylin Kozak have described such nucleotide sequences in 1984 (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC318541/). Kozak sequence can be different in different mRNAs and is was reported that different Kozak sequences influenced the translation level differently. In 2014, collective of scientists has published the data of direct  measurement of translation level for every possible Kozak sequence containing classic AUG start codon (https://pubmed.ncbi.nlm.nih.gov/25170020/) and computed the model of influence of the particular nucleotides on particular position in the Kozak sequence on the translation efficiency. 

It is known that human genome has a lot of variable positions and some of them are annotated as related with some diseases, some of them are referred as benign but the significance of some other genetic variants is uncertain for now. If the variant is located in the preotein-coding secuence or other well studied sequence, it is ok to predict the effect of such variant but if is is located in non-coding sequence the prediction becomes more unreliable. 

In this project I have tried to combine the data of Kozak sequence efficiency and genetic variants located in the Kozak sequences to predict possible pathogenicity of such variants to improve medical genetic analysis. 

My tasks were:
1. To select genetic variants from ClinVar and gnomAD databases, which are located in the Kozak sequences
2. To analyse the possible effect of these variants on the translation initiation efficiency
3. To train a model to predict the possible pathogenicity of the variants in the Kozak sequences
4. To create a tool for the correct annotation of such genetic variants

## Workflow description
### Requirements
 - python 3.8 
 - grep 3.4 
 - awk 5.0
 - bedtools 2.27 
 - R 4.1 
 - Jupyter notebook 6.0
### Data sources
 - gencode genome annotation and assembly - GRCh38.p13 Release 40 (https://www.gencodegenes.org/human/)
 - gnomAD variants with allele frequency >5% - kindly provided by Y.Barbitov
 - Clinvar variants - kindly provided by Y.Barbitov
 - Kozak sequence efficiency data - 10.15252/msb.20145136, supplementary
### Data preprocessing
#### Kozak sequence coordinates extraction
1. Filter with _grep_ and _awk_ the genome annotation to get only the protein-coding transcript and associated exons and CDS:

`zcat gencode.v40.annotation.gff3.gz | grep -v '^#' | awk -v FS="\t" -v OFS="\t" '$1 !~ "chrM"' |  awk -v FS="\t" -v OFS="\t" '$3 ~ "exon|transcript|CDS"' | awk -v FS=";" -v OFS=";" '$7 ~ "transcript_type=protein_coding"' > gencode.v40.annotation_withCDS.gff3 `

2. Run the script _Kozak_coordinates_finder.py_ with input file path and output file path as positional arguments

`python3 ./scripts/kozak_coordinates_finder.py gencode.v40.annotation_withCDS.gff3 gencode.v40_kozak_intervals_withCDSandstrands.bed` 

There are 4 columns in the resulted .bed file: 
 - chromosome
 - start position (0 based)
 - end position (0 based, it should be not included in the interval)
 - chain (+ or -)

#### Kozak sequence extraction from genome
1. Use _bedtools getfasta_ to extract the sequences according to the coordinates of the Kozak sequences

`bedtools getfasta -fi GRCh38.p13.genome.fa -bed gencode.v40_kozak_intervals_withCDSandstrands.bed > gencode.v40_kozak_seq.fasta`

In the resulted .fasta file were sequences of 12 nt (bedtools INCLUDES the right border of the interval). This file is usefull for us too, but we also need "pure" sequences (their lengths are 11 nt). 

2. Run the script _Kozak+1_to_Kozak.py_ which extract only 11 nt Kozak sequences and transform them to the reverse complement if they are located on '-' chain

`python3 ./scripts/Kozak+1_to_Kozak.py gencode.v40_kozak_seq.fasta gencode.v40_kozak_seq_pure.fasta` (!!!!! check!!!!!!)

#### Variants preparation
1. If the chromosome name in the .vcf is not 'chr1' but just '1', add 'chr' (there are 'chr1' in gencode and gnomAD data and '1' in clinvar and we need to unify them before intersecting)

`cat <filename>.vcf | awk '{print "chr"$0}' > <filename>_1.vcf`

2. Filter only the SNP (1 nt change to another 1 nt) (and I save the header too here)

`cat <filename>.vcf | grep '##' > <filename>_snp_with_header.vcf`

`cat <filename>.vcf | grep -v '^#' | awk -v FS="\t" -v OFS="\t" 'length($4) == 1'| awk -v FS="\t" -v OFS="\t" 'length($5) == 1' >> <filename>_snp_with_header.vcf`

These steps are similar for ClinVar and gnomAD .vcf files, so I use <filename> here to generalize (substitute it with your actual filename). 


### Main pipeline
1. Intersect the .bed Kozak coordinates with prepared .vcf files with _bedtools intersect_:

`bedtools intersect -a <filename>_snp_with_header.vcf -b gencode.v40_kozak_intervals_withCDSandstrands.bed -wa -wb > <filename>_in_gencode_kozak_CDSstrands_snps.txt`

There can be problems with intersection if the files are not sorted or chromosomes naming is not similar. 

The resulting files contains the full information about the variant located in the Kozak sequence (from .vcf) and about this Kozak sequence (from .bed). 
Unfortunately _bedtools intersect_ INCLUDES the right border of the interval too (as _bedtools getfasta_ in the **Preprocessing** section), so we have to filter the mutations which are not in the Kozak sequence actually. 

2. Filter the variants which position is not equal to the right interval border. 

`cat <filename>_in_gencode_kozak_CDSstrands_snps.txt | awk -v FS="\t" -v OFS="\t" '$2 != $11' > <filename>_in_gencode_kozak_CDSstrands_snps_noborders.txt`

bla blah

### Analysis with R
The 


### Model training



## Results




## Future plans
1. Fix the bugs:
 - To check ‘not AUG’ Kozak sequences (are they reliable or just script mistakes?) 
 - To analyse calculations for ClinVar variants more correctly
 - To compare the ClinVar variants and gnomAD variants
 - To train various models with corrected data (if it is possible)
2. Think
3. … and do some more science :)
