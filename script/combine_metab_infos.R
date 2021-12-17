library(ocean)
library(readr)
library(GSEABase)
library(metaboliteIDmapping)

metabolitesMapping <- as.data.frame(metabolitesMapping)

load("support/HMDB_metab_class.RData")
load("support/mapping_table.RData")

source("script/support_functions.R")

names(mapping_table)[1] <- "metabs"

metabs_df[grepl("HMDB[0-9]{9}",metabs_df$metabs),"metabs"] <- gsub("HMDB00","HMDB",metabs_df[grepl("HMDB[0-9]{9}",metabs_df$metabs),"metabs"])

full_metab_metadata <- merge(mapping_table, metabs_df, by = "metabs", all = T)

metabolitesMapping <- metabolitesMapping[which(metabolitesMapping$HMDB %in% full_metab_metadata$metabs),]
metabolitesMapping <- unique(metabolitesMapping[,c(5,8)])

# write(as.character(unique(full_metab_metadata$metabs)), file = "support/HMDB_IDs.txt")
names(metabolitesMapping) <- c("CID","metabs")

full_metab_metadata <- merge(full_metab_metadata,metabolitesMapping, by = "metabs")

write_csv(full_metab_metadata, file = "support/full_metab_metadata.csv")

recon3D_bigg_no_cofactor <- as.data.frame(
  read_csv("result/recon3D_bigg_no_cofactor.csv"))

recon3D_bigg_no_cofactor$source <- gsub("Metab__","cpd:",recon3D_bigg_no_cofactor$source)
recon3D_bigg_no_cofactor$target<- gsub("Metab__","cpd:",recon3D_bigg_no_cofactor$target)

recon3D_bigg_no_cofactor$source <- gsub("Gene[0-9]+__","",recon3D_bigg_no_cofactor$source)
recon3D_bigg_no_cofactor$target<- gsub("Gene[0-9]+__","",recon3D_bigg_no_cofactor$target)

metabolites <- rearrange_dataframe(recon3D_bigg_no_cofactor)

kegg_pathways <- import_gmt_pathways("support/c2.cp.kegg.v7.4.symbols.gmt")

metabolites_pathway_df <- map_pathways_to_metabolites(metabolites, pathways = kegg_pathways)
metabolites_pathway_df <- metabolites_pathway_df[grepl("HMDB",metabolites_pathway_df$metabolites),]
metabolites_pathway_df$metabolites <- gsub("_.*","",metabolites_pathway_df$metabolites)
metabolites_pathway_df$metabolites <- gsub("cpd:","",metabolites_pathway_df$metabolites)
metabolites_pathway_df <- unique(metabolites_pathway_df)

write_csv(metabolites_pathway_df, file = "result/HMDB_to_KEGG.csv")
write_csv(kegg_pathways, file = "result/genesymbol_to_KEGG.csv")
