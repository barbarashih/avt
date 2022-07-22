# Amplicon variance tracker (avt)

This tool tracks microbiome variance and identifies bacterial genera responsible for significant variance shifts using 16S rRNA gene metabarcoding for longitudinal data. The unique value of this approach is that it goes beyond presence, absence and abundance of genera, as is commonly the case in 16S microbiome studies. Since gut and lung microbiome analyses use DNA extracted from faecal and sputum samples respectively, they only represent a distribution of their respective anatomical compartments. So, making inferences based on presence and absence of individual units (genera) of the distribution might be misleading. In this prototype tool, variance is tracked using three microbiome attributes/parameters i.e. prevalence, abundance and co-occurrence network along a temporal scale identifying genera responsible for dramatic shifts in variance at each temporal point.

## Install

### Requirements:
<ul>
<li>R (> version 4.1.2)</li>
<li>devtools (R library)</li>
<li>phyloseq (R library)</li>
<li>mirobiome (R library)</li>
</ul>

#### Install pre-requisits
        install.packages("devtools")

        if (!require("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install("phyloseq")
        BiocManager::install("mirobiome")

#### Install avt from github
        install_github("barbarashih/avt")

### Check the installation with the test dataset
Download the test data folder ("avt_test_data") from this github.

Set the working directory to avt_test_data.

        setwd("avt_test_data")

Organise the avt inputs.

        organised_data <- avt_data(data_dir="SRP221070", grouping_column="Sampling.point")

Find the top genera. If your data has a large number of genera, you may want to consider setting centrality_p, centrality_r_upper and centrality_r_lower to filter Pearson correlation coefficients p and r values for the centrality calculation.

        top_genera <- topGenera(avt_dat=organised_data)

Calculate the change in variace when each genus is removed from the top genera. This step can take a little while to run. 

        top_genera_variance <- calculateVariance(avt_dat=organised_data, top_genera=top_genera)

Plot the change in variance.

        plotEigeSum(top_genera_variance) 

Print a summary on variance change with the removal of each genus.

        variance_change_summary <- orderedVarianceChange(top_genera_variance)

The default output folder name is "avt_out". 
