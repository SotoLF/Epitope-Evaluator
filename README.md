# Epitope-Evaluator
**Epitope Evaluator** is a Shiny app to analyze predicted T-cell epitopes. Epitope-Evaluator includes 6 tools: 
1. Epitope Distribution
2. Epitope Intersection
3. Epitope Density
4. Epitope Location 
5. Epitope Promiscuity
6. Epitope Conservation
 
## Input File
Epitope-Evaluator requires as input: 1) a multi-FASTA file containing the IDs and the sequences of the antigens, and 2) the prediction file obtained from a T-cell epitope predictor (Figure 1D). The user must indicate the predictor being used from a set of options or indicate “others” if not obtained from any of the listed predictors. In this case, the prediction file should have the following columns: the peptide sequence, its position within the protein, the protein ID, the protein length, and subsequent columns corresponding to each of the MHC alleles evaluated, where each value indicates a score for each epitope (Figure 1E).  In addition to this, users must indicate whether the score in the table corresponds to the “percentile rank” or “binding affinity score”. The **percentile rank** score is a number between 0 and 100. It is defined as the rank of the predicted binding score compared to a set of random natural peptides by several predictors of T-cell epitopes [(32711842](https://pubmed.ncbi.nlm.nih.gov/32711842/) , [32406916](https://pubmed.ncbi.nlm.nih.gov/32406916/), [32308001)](https://pubmed.ncbi.nlm.nih.gov/32308001/). The applicative will automatically identify whether the epitopes are class I or class II based on the name of the MHC alleles. 

The title section indicates the different tools that are available in the web tool. Each of these tools is independent, thus, users can run all analyses simultaneously. Each of the six tools and the ‘Input’ tab has four different sections: 
1. The parameters section : Located on the left side of the web application, allows users to set different options and parameters for the corresponding tool
2. The title section :
3. The help section : It describes the functionality of each tool, details each parameter, and explains the plots and tables returned in the output section.
4. The output section : It shows the plots and tables from the selected analyses which are downloadable. All the tools contain interactive plots where users can zoom in/out, select regions, and obtain more information by hovering

<p align="center">
 Representation of Epitope-Evaluator and its tools
</p>
<p align="center">
 <img src="https://github.com/SotoLF/Epitope-Evaluator/blob/main/Images/Github1.PNG">
</p>


Each of the tools shows 
1. Graphics that can be downloaded as png files by clicking on the camera icon **Download plot as png**
2. Tables that can be downloaded by clicking on the **Download table** button.

The R scripts of Epitope-Evaluator can be freely downloaded and launched locally.



# Epitope Distribution: Distribution of epitopes by MHC allele

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
<p align="center">
 Histogram and Results table of the Spike protein Epitope Distribution for MHCII 
</p>
<p align="center">
 <img width="400" height="500" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_Epitope_Distribution.png?raw=true">
 <img width="450" height="500" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_table_Epitope_Distribution.png?raw=true">
</p>

<p align="center">
 Duplicate Epitopes of the Spike protein Epitope Distribution for MHCII
</p>
<p align="center">
 <img src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_Duplicate_Epitopes.png?raw=true">
</p>

# Epitope Density: Correlation between number of epitopes and protein lengths

Epitope Density shows a scatter plot where the protein length in amino acids is on the x-axis and the number of epitopes is on the y-axis. Hovering over each point shows the name of the protein, number of epitopes, length of the protein, and the epitope density. The epitope density is calculated as the number of epitopes divided by the length of the proteins. It also shows a heatmap representing the number (or density) of epitopes in each protein across each allele. 
To obtain the scatter plot, users must indicate the maximum cutoff %rank to consider a peptide as an epitope and click on the **Run analysis** button. This will produce a scatter plot and a table containing the ID of the proteins, their length in amino acids, the number of epitopes, and the epitope density. Clicking on any row of the table will highlight the point in the scatter plot. 
To obtain the heatmap, in addition to the maximum cutoff%rank, users must select the plot-type (heatmap or bar plot) and the fill-type (by number or density of epitopes). We recommend using the plot-type **heatmap** when several proteins are analyzed. When heatmap is selected, this shows the alleles on the x-axis and the proteins on the y-axis. The color intensity indicates the log 10 of the number or the density of epitopes. When **barplot** is selected, each protein is represented as a bar while alleles are on the x-axis and the number (or density) of epitopes are on the y-axis. Hovering over each bar (or cell in the heatmap) will show the protein ID, the allele, the number of epitopes, the length of the protein, and the density of epitopes. 

PARAMETERS
- **Cutoff % rank:** Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10).
- **Color by:** Fill heatmap by number or density of epitopes.(Default = Epitopes Number)
- **Plot type:** Heatmap is recommended for several proteins, while bar plot is recommended for few proteins. (Default = Heatmap)

RESULTS
<p align="center">
 Plot and Results table of the Spike protein Epitope Density for MHCII 
</p>
<p align="center">
 <img width="400" height="400" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_Epitope_Density.png?raw=true">
 <img width="500" height="300" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_table_Epitope_Density.png?raw=true">
</p>

<p align="center">
 Heatmap of the Spike protein Epitope Density for MHCII
</p>
<p align="center">
 <img width="800" height="350" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_heatmap_Epitope_Density.png?raw=true">
</p>
<p align="center">
 Barplot of the Spike protein Epitope Density for MHCII
</p>
<p align="center">
 <img width="800" height="400" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_barplot_Epitope_Density.png?raw=true">
</p>

# Epitope Location: Amino acid position of the epitopes in proteins

Epitope Location shows the location of epitopes within the protein (blue bar). Epitopes are represented on a color gradient from yellow to red reflecting the number of MHC alleles they are predicted to bind to. Hovering on each epitope will show the amino acid sequence, the position within the protein, and the MHC alleles that are predicted to recognize it. Users must select the protein, set a cutoff %rank to filter epitopes, and indicate the MHC alleles of interest. When multiple alleles are selected, users also need to specify if they are interested in epitopes that are predicted to bind all the selected alleles (AND) or at least one of the selected alleles (OR). This tool also shows a table containing the sequence of the epitope, the start and end amino acid position within the protein, the alleles, and the number of alleles to which they were predicted to bind. This analysis may take a few minutes depending on the number of epitopes and the length of the protein. 

PARAMETERS
- **Protein:** List of proteins obtained from the input file. Users must select a protein to be analyzed. (Default = first protein)
- **Cutoff % rank:** Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10).
- **Condition over MHC alleles:** When multiple alleles are selected, users must indicate if they are interested in shared epitopes that bind all the selected MHC alleles (AND) or epitopes that bind at least one of the MHC alleles selected (OR).
- **MHC alleles:** List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = first allele)

RESULTS
<p align="center">
 Plot and Results table of the Spike protein Epitope Location for MHCII
</p>
<p align="center">
 <img width="650" height="250" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_Epitope_Location.png?raw=true">
 <img width="300" height="300" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_table_Epitope_Location.png?raw=true">
</p>

# Epitope Promiscuity: Number of alleles and binding rank per epitopes

Epitope Promiscuity shows epitopes predicted to bind a minimum number of MHC alleles. The tool depicts a heatmap where alleles are on the x-axis and epitopes are on the y-axis. Users should set the strong and weak binding cutoff %rank to filter the strong binder and weak binder epitopes. Strong binder epitopes (SB) are the epitopes with a %rank less or equal than the strong binding cutoff %rank. Weak binder epitopes (WB) are the epitopes with a %rank greater than the strong binding cutoff %rank but less or equal than the weak binding cutoff %rank. In the heatmap, SB and WB epitopes are colored in red and orange, respectively. Users also need to select a minimum number of MHC alleles to filter and show only epitopes that are predicted to bind to at least this number of alleles. This tool also shows a table containing the amino acid sequence, the position within the protein, the ID of the protein, and the number of alleles to which they are predicted to bind. 

PARAMETERS
- **Minimum number of alleles:** This indicates the minimum number of MHC alleles that an epitope needs to bind to be shown in the heatmap.
- **Strong Binding Cutoff % rank:** Epitopes with a %rank lower or equal than this cutoff are considered as strong binder epitopes (SB).
- **Weak Binding Cutoff %rank:** Epitopes with a %rank lower or equal than this cutoff but greater than Strong Binding Cutoff are considered weak binder epitopes (WB).

RESULTS
<p align="center">
 Plot and Results table of the Spike protein Epitope Promiscuity for MHCII
</p>
<p align="center">
 <img width="400" height="500" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_Epitope_Promiscuity.png?raw=true">
 <img width="400" height="500" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_table_Epitope_Promiscuity.png?raw=true">
</p>


# Epitope Intersection: Set of epitopes shared between different MHC allele combinations

Epitope Intersection shows an Up-Set plot where the MHC alleles are represented as the sets and the epitopes are represented as the elements. The selected MHC alleles are represented with bars on the left side indicating the number of predicted epitopes for each MHC allele. The number of epitopes for each MHC allele or intersection of MHC alleles is represented with bars at the top. Individual points in the grid indicate epitopes binding to a specific MHC allele while connected points indicate epitopes that can bind to multiple MHC alleles. This tool also provides a table containing the number of alleles, the alleles, the epitope sequences within each region, and their respective number of epitopes. 

PARAMETERS
- **Cutoff %rank:** Maximum %rank to consider a peptide as predicted epitope. (Default MHC Class I = 2, MHC Class II = 10)
- **MHC alleles:** List of MHC alleles obtained from the input file. Users can choose more than one allele. (Default = first 5 alleles)

RESULTS
<p align="center">
 Plot of the Spike protein Epitope Intersection for MHCII
</p>
<p align="center">
 <img width="800" height="500" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_Epitope_Intersection.png?raw=true">
</p>
<p align="center">
 Results table of the Spike protein Epitope Intersection for MHCII
</p>
<p align="center">
 <img width="900" height="600" src="https://github.com/SotoLF/Epitope-Analyzer/blob/main/Examples/Results/Spike_MHCII_table_Epitope_Intersection.png?raw=true">
</p>

# How to cite Epitope Analyzer
Paper Citation

# Contact
Juan Fuxman Bass - fuxman@bu.edu 

Luis F Soto - lufesu98@gmail.com
