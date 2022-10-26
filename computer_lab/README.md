# Computer lab

This folder contains data and R scripts to be used in the computer lab.
The data is 3 years from the Isa time series (2011-2014) processed through USearch with clustering at 99%,
arbitrarily removing any OTU with less than 50 reads in the total dataset and removing OTUs assigned to Streptophyta, Metazoa and Fungi.
Rarefied to 5398 reads per sample. Taxonomic assignment with PR2 using blast. Only DNA from 25m is included.

There are three main script
- [Computer lab I](Computer_lab_I.Rmd)
- [Computer lab II](Computer_lab_II.Rmd)
- [Computer lab III](Computer_lab_III.Rmd)

Example scripts:
[Barplot_examples.R](Barplot_examples.R) shows how you can subset data and make barplots of the Otu table. 
