# Epitope-Analyzer
Epitope Analyzer is a Shiny app aimed to filter and analyze predicted T-cell epitopes. Epitope Analyzer includes 5 tools: 
1. Epitope Distribution
2. Epitope Density
3. Epitope Location
4. Epitope Promiscuity 
5. Epitope Intersection
 
It needs a txt format file tab-delimited where the first 4 columns are the peptide sequence, the position of the peptide within the protein, the length of the protein, and the name of the protein. The following columns should be named with the allele containing the %rank score assigned to each peptide for each MHC allele. The %rank score is a number between 0 and 100. It is defined as the rank of the predicted binding score compared to a set of random natural peptides by several predictors of T-cell epitopes ( 32711842 32406916 32308001 ). 

Usually, the output of these predictors requires to be parsed before uploading in the shiny app. Once the input file is uploaded, parameters need to be set by the user followed by clicking on the Run analysis or Plot buttons. 
Each of the tools shows 
1. Graphics that can be downloaded as png files by clicking on the camera icon ( Download plot as png)
2. Tables that can be downloaded by clicking on the Download table button.

We have included 2 example files that can be downloaded and used as input files on Epitope Analyzer.

- File 1 : Predicted MHC Class I epitopes from SARS-CoV-2.
- File 2 : Predicted MHC Class II epitopes from SARS-CoV-2

The R scripts of Epitope Analyzer can be freely downloaded from GitHub and launched locally.

# How to cite Epitope Analyzer
Paper Citation

# Contact
Juan Fuxman Bass - fuxman@bu.edu 

Luis F Soto - lufesu98@gmail.com
