library(ocean)
library(decoupleR)
library(readr)

HMDB_to_KEGG <- as.data.frame(
  read_csv("result/HMDB_to_KEGG.csv"))
genesymbol_to_KEGG <- as.data.frame(
  read_csv("result/genesymbol_to_KEGG.csv"))



network <- metabolites_pathway_df    #input 1: pathways and their associated metabolites
network$mor <- 1                     #add new column for mode of regulation
network$likelihood <- 1              #add new column for edge likelihood


#mat should be a dataframe with features (metabolite or genes) as row names, 
#and each column being the correspoding stat of a given sample or contrast.
#the feature names should be the same as in network

##Perform enrichment analysis
enrichment <- run_wmean(mat, network,
                        .source = .data$pathway,
                        .target = .data$metabolites,
                        .mor = .data$mor, 
                        .likelihood = .data$likelihood,  
                        times = 1000)  #randomize matrix