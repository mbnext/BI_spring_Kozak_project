# Prediction of pathogenicity of genetic variants in Kozak sequences

Student: Marianna Baranovskaia

Supervisors: Yury Barbitov (Bioinformatics Institute), Michail Skoblov (Research Center of Medical Genetics)

## Introduction
Kozak sequence is a consensus nucleotide environment of the start codon in the most of the eukaryotic mRNAs, involved in the translation initiation. Marylin Kozak have described such nucleotide sequences in 1984 [1] . Kozak sequence can be different in different mRNAs and is was reported that different Kozak sequences influenced the translation level differently. In 2014, collective of scientists has published the data of direct  measurement of translation level for every possible Kozak sequence containing classic AUG start codon and computed the model of influence of the particular nucleotides on particular position in the Kozak sequence on the translation efficiency [2]. 

It is known that human genome has a lot of variable positions and some of them are annotated as related with some diseases, some of them are referred as benign but the significance of some other genetic variants is uncertain for now. If the variant is located in the protein-coding secuence or other well studied sequence, it is ok to predict the effect of such variant but if is is located in non-coding sequence the prediction becomes more unreliable. 

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
 - gencode genome annotation and assembly - GRCh38.p13 Release 40 [3]
 - Clinvar variants - extracted data from [4] were kindly provided by Y.Barbitov
 - gnomAD variants with allele frequency >5% - extracted data from [5] were kindly provided by Y.Barbitov 
 - Kozak sequence efficiency data - [2], supplementary
### Data preprocessing
#### Small useful functions
I have write small code in the script _dna_rna_functions.py_:
1. _dna_rev_com(seq)_ to transform DNA into reverse complement DNA
2. _dna_to_rna_convert(seq)_ to transform DNA into RNA
3. _mutation_type(ref_triplet, alt_triplet)_ to define the type of mutation (synonymous, missense or nonsense). 
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

2. Run the script _Kozak+1_to_Kozak.py_ which extract only 11 nt Kozak sequences and transform them to the reverse complement if they are located on '-' chain. There are 3 positional arguments here: output file path, .fasta with extracted Kozak sequences (containing the additional 1 nt), extracted Kozak coordinates in .bed file. 

`python3 ./scripts/Kozak+1_to_Kozak.py gencode.v40_kozak_seq_pure.fasta gencode.v40_kozak_seq.fasta gencode.v40_kozak_intervals_withCDSandstrands.bed` 

#### Variants preparation
1. For gnomAD file: filter only the SNP (1 nt change to another 1 nt) (and I save the header too here)

`cat <filename>.vcf | grep '##' > <filename>_snp_with_header.vcf`

`cat <filename>.vcf | grep -v '^#' | awk -v FS="\t" -v OFS="\t" 'length($4) == 1'| awk -v FS="\t" -v OFS="\t" 'length($5) == 1' >> <filename>_snp_with_header.vcf`

2. For ClinVar file add 'chr' to the chromosome naming if it has no 'chr' and perform the same filtration.

`cat <filename>.vcf | grep '##' > <filename>_snp_with_header.vcf`

`cat <filename>.vcf | grep -v '^#' | awk '{print "chr"$0}' | awk -v FS="\t" -v OFS="\t" 'length($4) == 1'| awk -v FS="\t" -v OFS="\t" 'length($5) == 1' >> <filename>_snp_with_header.vcf`


### Main pipeline
1. Intersect the .bed Kozak coordinates with prepared .vcf files with _bedtools intersect_:

`bedtools intersect -a <filename>_snp_with_header.vcf -b gencode.v40_kozak_intervals_withCDSandstrands.bed -wa -wb > <filename>_in_gencode_kozak_CDSstrands_snps.txt`

There can be problems with intersection if the files are not sorted or chromosomes naming is not similar. 

The resulting files contains the full information about the variant located in the Kozak sequence (from .vcf) and about this Kozak sequence (from .bed). 
Unfortunately _bedtools intersect_ INCLUDES the right border of the interval too (as _bedtools getfasta_ in the **Preprocessing** section), so we have to filter the mutations which are not in the Kozak sequence actually. 

2. Filter the variants which position is not equal to the right interval border. 

`cat <filename>_in_gencode_kozak_CDSstrands_snps.txt | awk -v FS="\t" -v OFS="\t" '$2 != $11' > <filename>_in_gencode_kozak_CDSstrands_snps_noborders.txt`

3. Calculate the position in Kozak where the substitution occurs. I use here numbers from 0 to 10 to describe the positions: 

classic notation || ordered position (AUG codon is dedicated)
  
-6 || 0
  
-5 || 1
  
-4 || 2
  
-3 || 3
  
-2 || 4
  
-1 || 5
  
+1 || 6 A
  
+2 || 7 U
  
+3 || 8 G
  
+4 || 9
  
+5 || 10
  
I decided to add the column according to the following rules: 
  - if "+" in $12, $13 = $2-$10-1 (variant position minus Kozak start position minus 1 because .bed is 0 based)
  - if "-" in $12, $13 = $11-1-($2+1)-1 = $11-$2-1 (Kozak end position minus 1 (because 'Kozak end' here is actually Kozak start on the '-' chain and right border is not included) minus variant position and 1 (because .bed is 0 based) and minus 1 again because my notation is 0 based)

Additionally I delete here the double chromosome column ($9). 

`cat <filename>_in_gencode_kozak_CDSstrands_snps_noborders.txt | awk -v FS="\t" -v OFS="\t" '$12 ~ "+"' |  awk -v FS="\t" -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $10, $11, $12, $13 = $2-$10-1}' > <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos.txt`
  
`cat <filename>_in_gencode_kozak_CDSstrands_snps_noborders.txt | awk -v FS="\t" -v OFS="\t" '$12 ~ "-"' |  awk -v FS="\t" -v OFS="\t" '{print $1, $2, $3, $4, $5, $6, $7, $8, $10, $11, $12, $13 = $11-$2-1}' >> <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos.txt`
  
`sort -k1,1 -k2,2n <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos.txt > <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted.txt`

4. Run the script for 'crude' annotation _mutation_sense_finder.py_: it analyses the variant position within the Kozak sequence and classifies the variant to the one of 5 classes: 'upstream' (the position is in range from 0 to 5), 'no_start' (the position is in range from 6 to 8, in the start codon), 'synonymous'/'missense'/'nonsense' (the position is 9 or 10 and Kozak sequences plus 1 nt are used to understand the change in 2nd codon). Before the annotation the script checks if the Ref letter from the .vcf file is on the dedicated position in the Kozak sequence extracted from the genome (if FALSE, the variant gets the class 'Error in annotation'). There are 4 positional arguments in the script: input file path, output file path, .fasta with extracted Kozak sequences (+1 nt) and the path to a file where the script writes data about 'Error in annotation'. 
  
`python3 ./scripts/mutation_sense_finder.py <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted.txt <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted_annot.txt gencode.v40_kozak_seq.fasta gencode.v40_kozak_intervals_withCDSandstrands.bed  mutation_sense_errors.txt`
  
The scipt add also the column with the type of Kozak sequence related with the variant: 'AUG_Kozak' or 'not_AUG_Kozak' (if there is no 'AUG' codon in position 6-8). 'not_AUG_Kozak' sequences are not involved in the analysis because they do not match to the Kozak sequences listed in the paper of Noderer at al. 

5. Run the script for accurate annotation of upstream and synonymous mutation _upstream_mutation_annotator.py_. The script reconstitutes the reference Kozak sequence and finds its efficiency in the list of Kozak sequences from the Noderer's paper and its confidence interval, then reconstitutes the alternative Kozak sequence and finds its efficiency and its confidence interval, calculated the relative efficiency (ratio of Alt. Kozak efficiency to Ref. Kozak efficiency) and add to the file the corresponding columns. There are 4 positional arguments in the script: input file path, output file path, .fasta of 11 nt Kozak sequences extracted from the genome and .txt file with all Kozak sequences efficiencies. 
  
`python3 ./script/upstream_mutation_annotator.py <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted_annot.txt <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted_annot_combo.txt gencode.v40_kozak_seq_pure.fasta Kozak_efficiency_rawdata.txt`
  
If the variant has an annotation 'no_start', 'Error in annotation', 'missense'/'nonsense' or Kozak type 'not_AUG_Kozak', there are '.' in the added columns. 
  
6. If we analyse the ClinVar file, we can extract the information about clinical significance and gene name from the column '##INFO' with the script _clnsig_finder.py_:
  
`python3 ./scripts/clnsig_finder.py <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted_annot_combo.txt clinvar_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted_annot_combo_clnsig.txt`
  
There are 2 positional arguments in the script: input and output files pathways. It adds 2 columns to the file: cln_sig and gene.
  
  In provided gnomAD file was no information about clinical significance and gene name, so I have added the columns 'cln_sig' with the value 'Benign' (probably it was a mistake) and 'gene' with '.' (to maintain the column number in clinvar and gnomAD file). 
  
7. Let's make some corrections: 

`cat <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted_annot_combo_clnsig_1.txt  | awk -v FS="\t" -v OFS="\t" '{print $3, $1, $2, $4, $5, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24}' > <filename>_in_gencode_kozak_CDSstrands_snps_noborders_pos_sorted_annot_combo_clnsig_dataset.txt`
  
  and add the header of the dataset manually:
"ID	chromosome	position	Ref	Alt	Kozak_start	Kozak_end	Chain	Kozak_variant_position	variant_annotation	Kozak_type	Ref_Kozak_efficiency	Ref_Kozak_lower	Ref_Kozak_upper	Alt_Kozak_efficiency	Alt_Kozak_lower	Alt_Kozak_upper	Change_description	Relative_efficiency	Clin_Sig	Gene"

### Analysis with R
The datasets for ClinVar and gnomAD variants were downloaded as data frames and then joined with _rbind()_. Vizualization was performed with the package _ggplot2_. The _reshape2_ package has to be downloaded too. The details are in the file **spring_project.Rmd** (knitted version is [./spring_project.pdf](https://github.com/mbnext/BI_spring_Kozak_project/blob/draft/spring_project.pdf))

### Model training
The dataset for ClinVar variants was downloaded in Jupyter notebook. I have tried to learn the _DecisionTreeClassifier_ from the library _Sci-kit Learn_ but the results were diffucult to interpret and they are not represented here. May be we will return to this task after corrections in the previous described analysis. 

## Results
Total number of the sequences and variants during the pipeline:
 - Extracted Kozak sequences from genome assembly: ~39000
 - Total gnomAD and Clinvar variants intersected with extracted Kozak’s: 7921 
 - Variants which can be referred as benign or pathogenic and likely pathogenic (variants from gnomAD with AF>5% were referred as benign here): 2657
 - Variants located in Kozak sequences with classic AUG start codon: 1984
 - Variants with ‘upstream’ and ‘synonymous’ annotation: 807
 - Upstream and synonymous variants which have nonintersected confidence intervals between Ref and Alt sequences: 92

We can see that only a few variants have significant effect of the mutation. The suspicious facts are that a lot of 'not_AUG_Kozak' sequences were in the extracted from the genome Kozak set, and the second fact is that about 500 variants from gnomAD was intersected (from about 8000000 of SNP), it seems to be too small number (maybe there is a mistake somewhere). 


![plot_F_04](https://user-images.githubusercontent.com/34398217/169703893-dd993030-71f5-4a46-8b79-37c95367b439.jpeg)

**Fig. 1. Distribution of the annotation types among the 2657 variants of known pathogenicity (AUG + not_AUG Kozak's)**


![plot_F_09](https://user-images.githubusercontent.com/34398217/169703902-68348e85-e821-40c9-968b-f5659cd2ccdf.jpeg)

**Fig. 2. Distribution of the annotation types among the 1984 variants of known pathogenicity in AUG Kozak's**


![plot_F_09_notAUG](https://user-images.githubusercontent.com/34398217/169713477-729d7b97-0cf0-4637-9dda-4b81da426b42.jpeg)

**Fig. 3. Distribution of the annotation types among the 673 variants of known pathogenicity in not-AUG Kozak's**

Distributions in AUG and not-AUG Kozak sequences are different. 

![plot_F_19](https://user-images.githubusercontent.com/34398217/169703919-a5cf8add-64ce-4980-8c33-c8f7369bfb8a.jpeg)

**Fig. 4. Ref. and Alt. Kozak sequences efficiencies for pathogenic and benign variants (for 92 ‘significant’ variants)**

![plot_F_22](https://user-images.githubusercontent.com/34398217/169713391-8f8b7ff5-1b4c-4aa4-aa6b-b15f6274664d.jpeg)

**Fig. 5. Distribution of relative efficiencies (Alt.Eff/Ref.Eff.) of Kozak sequences  for pathogenic and benign variants according to the variant position in upstream Kozak sequence (for 92 ‘significant’ variants)** 
One outlier (near 4) was deleted from the plot. 


More plots are in the file [./spring_project.pdf](https://github.com/mbnext/BI_spring_Kozak_project/blob/draft/spring_project.pdf).

## Conclusions
According to performed analysis, there are no significant differences in the predicted translation efficiency of the reference and alternative Kozak sequences referred as benign or pathogenic (but the particular position can be significant, the additional more carefull analysis seems to be performed after check of all of the scripts and calculation procedures).

## Future plans
1. Fix the bugs:
 - To check ‘not AUG’ Kozak sequences (are they reliable or just script mistakes?) 
 - To analyse calculations for ClinVar variants more correctly
 - To compare the ClinVar variants and gnomAD variants
 - To train various models with corrected data (if it is possible)
2. Think
3. … and do some more science :)

## References
1. Kozak M. Compilation and analysis of sequences upstream from the translational start site in eukaryotic mRNAs. Nucleic Acids Res. 1984;12(2):857-872. doi:10.1093/nar/12.2.857 
2. Noderer WL, Flockhart RJ, Bhaduri A, et al. Quantitative analysis of mammalian translation initiation sites by FACS-seq. Mol Syst Biol. 2014;10(8):748. Published 2014 Aug 28. doi:10.15252/msb.20145136
3. GRCh38.p13 Release 40, URL: https://www.gencodegenes.org/human/
4. ClinVar database, URL: https://www.ncbi.nlm.nih.gov/clinvar/
5. gnomAD database, URL:https://gnomad.broadinstitute.org/
