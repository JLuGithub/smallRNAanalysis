# smallRNAanalysis


# About our Scrips 

Innate and exogenous dsRNAs have been reported to be targets for mammalian Dicer, and thus be cleaved into small-interfering RNAs (siRNAs). Here we proveded this script for siRNA analysis, and can be applied on viral siRNA (vsiRNAs), transfected dsRNAs, or exogenous dsRNAs. The script would analyze the patten of paired small RNAs and report several metrics that help users confirm the generation of siRNAs from desired templates.
Our script has the following features:

-  One step generation of all results with required input:
     - User friendly operations, especially for wetters. 
     - All results of same normalization method stored in one folder.
     - Intermediated data stored in the xls files and can be loaded for further analysis. 
- The script would generate a compressed intermediate fastq for storage or transfer, but capable for re-performing all analysis.
- The input files are easy to prepare.

--------------------------

# Installation

We have provided yaml file and can be easily installed by conda or mamba.

`conda env create -f small_RNA.yml -n `

**Note:** MACOS 11 or higher is pre-installed with perl so the yaml files doesn't contain this installing information, other users may need to install perl manually. 

You can also install all softwares by the following steps.
```{bash}
conda env create -n exo
sudo apt-get install perl
sudo apt-get install libgd-dev
sudo perl –MCPAN –e shell
install GD
install Spreadsheet::WriteExcel
sudo apt-get install bowtie
sudo apt-get install libtbb-dev
```


# Input file preparation
**Please prepare the following files in the working directory:**
- Fastq files of small RNA library (sequecing on SE50 mode is preferred,if the sequencing is performed on PE150 mode, please use R1.fastq as input file) `${YourRawReads}.fastq`
- Reference file for miRNAs in .fa format, naming as `${ref_miRNAs}.miRNA.fa `
- Reference file for exogenous dsRNAs (or viral genome) in .fa format, naming as `${ref_dsRNA}.exo.fa `
- The script `Autoanalysis.perl`

**If you don't have the miRNA file of your species, please generate reference file by:**

Download mature miRNA reference of all species from miRbase(https://www.mirbase.org/download/mature.fa)
Use the following commands to extract the miRNA information of species that you need
```{bash}
All_mature_file = mature.fa
specie = pal 
```
In our study, we used ‘pal’ for ***P**teropus **a**lecto*, ‘mmu’ for ***M**us **mu**sculus*, and ‘hsa’ for ***H**omo **sa**piens*, the naming convention can be applied to other species.. 
If the desired species are absent in the file, please use relatived species instaed.
```
awk -v sp="$species" '/^>/{p=($0 ~ sp)} p' $ All_mature_file > ${species}.miRNA.fa
```
Now you can see the miRNA reference in the new generated file
```
ls ${species}.miRNA.fa
# >pal-miR-1-5p
# ACATACTTCTTTATGTACCCATA
# >pal-miR-1-3p
# TGGAATGTAAAGAAGTATGTAT
# >pal-miR-1-2-5p
# ACATACTTCTTTATATGCCCAT
```



# Usage
**Please prepare the following files in the working folder:**
- Go to your working folder:
`cd ${working_dir}`
- Activate the enviroment
`conda activate exo`
- Run the script
`perl Autoanalysis.perl`
- Wait for few minutes, you will see all of these reports:
```
# Testdata.raw.fastq is the input fastq file
# The file is not trimmed, autotrimming running
# 5' none and 3' AGATCGGA is default setting, see line 2388
# Trimming, please wait
# The mapping will be done allowing no mismatches.

# dsRNA.exo.fa is the input reference file of exogenouse dsRNA
# Segment 1 : >mydsRNA, complete sequence 257 nt
# maturepal22.miRNA.fa is the input reference file of miRNA 

# processing no normed data
# Done.
# processing data normalized by total 18_28nt
# Done.
# processing data normalized by miRNA
# Done.

# Program done. Please find all results in 'No_Norm' 'Norm_to_total' and 'Norm_to_miRNA' folder.
```

# Cite

**RNA Interference Functions as an Antiviral Immunity Mechanism in Mammals**
(Algorithm for exogenous small-interfering RNA analysis) 
`https://www.science.org/doi/10.1126/science.1241911`

**Protocol for identifying Dicer as dsRNA binding and cleaving reagent in response to transfected dsRNA**
(Experiments and analysis for study of exogenous dsRNA generated small-interfering RNA)
`Link to STAR protocol`


# Contact

Please contact yli@ioz.ac.cn or jl5103@cumc.columbia.edu or raise an issue in the github repo with any questions about installation or usage. 
