# Epitope-Analyzer
**Epitope Analyzer** is a Shiny app aimed to filter and analyze predicted T-cell epitopes. Epitope Analyzer includes 5 tools: 
1. Epitope Distribution
2. Epitope Density
3. Epitope Location
4. Epitope Promiscuity 
5. Epitope Intersection
 
It needs a txt format file tab-delimited where the first 4 columns are the peptide sequence, the position of the peptide within the protein, the length of the protein, and the name of the protein. The following columns should be named with the allele containing the %rank score assigned to each peptide for each MHC allele. The **%rank** score is a number between 0 and 100. It is defined as the rank of the predicted binding score compared to a set of random natural peptides by several predictors of T-cell epitopes [(32711842](https://pubmed.ncbi.nlm.nih.gov/32711842/) , [32406916](https://pubmed.ncbi.nlm.nih.gov/32406916/), [32308001)](https://pubmed.ncbi.nlm.nih.gov/32308001/). 

Usually, the output of these predictors requires to be parsed before uploading in the shiny app. Once the input file is uploaded, parameters need to be set by the user followed by clicking on the **Run analysis** or **Plot** buttons. 
Each of the tools shows 
1. Graphics that can be downloaded as png files by clicking on the camera icon **Download plot as png**
2. Tables that can be downloaded by clicking on the **Download table** button.

We have included examples files that can be downloaded [here](https://github.com/SotoLF/Epitope-Analyzer/tree/main/Examples) and used as input files on Epitope Analyzer.

The R scripts of Epitope Analyzer can be freely downloaded and launched locally.

# Distribution of epitopes by MHC allele

Epitope Distribution shows:
1. The distribution of the number of unique epitopes that were predicted to bind to selected MHC alleles within an %rank interval 
2. A heatmap showing the presence of epitopes across different proteins. 

To obtain the distribution of epitopes, users should select one or more MHC alleles, set the min %rank and the max%rank to consider a peptide as an epitope, and select the plot type which could be a histogram or cumulative histogram. The histogram and cumulative histogram plots show the %rank intervals on the x-axis and the number of epitopes on the y-axis. The density is represented with a red line in the histogram plot. To show the epitopes present across multiple proteins, users need to indicate a minimum number of proteins in which the epitopes are contained. The Shiny app will show a heatmap where the epitopes are on the x-axis while proteins are on the y axis. The cells are colored with blue or white indicating the presence or absence of the epitope, respectively. 

PARAMETERS
- **Relation among MHC alleles:** When multiple alleles are selected, users must indicate if they are interested in shared epitopes that bind all the selected MHC alleles (AND) or epitopes that bind at least one of the MHC alleles selected (OR).
- **MHC alleles:** List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = first allele)
- **Min % rank:** Minimum %rank to filter the epitopes. (Default = 0)
- **Max % rank:** Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10).
- **Step:** Width of the bins in the histogram plot. (Default = 0.1)
- **Plot type:** The distribution of epitopes can be shown as a histogram or a Cumulative histogram. If histogram is selected, the density is represented with a red line.
- **Min number of proteins:** Minimum number of proteins where epitopes are conserved (Default = 2).

RESULTS


![](https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png){width='100px'}



# How to cite Epitope Analyzer
Paper Citation

# Contact
Juan Fuxman Bass - fuxman@bu.edu 

Luis F Soto - lufesu98@gmail.com
