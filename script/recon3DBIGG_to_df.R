library(R.matlab)
library(stringr)
library(readr)
library(ocean)

recon3D_BIGG <- readMat("Recon3D.mat")


recon3D_paper <- readMat("Recon3DModel_301.mat")

recon3D_BIGG <- recon3D_BIGG$Recon3D
recon3D_paper <- recon3D_paper$Recon3DModel

#get the stochio matrix
S <- as.matrix(recon3D_BIGG[[11]])

#get the list of gene rules for reactions
reaction_list <- recon3D_BIGG[[7]]

#get reversible reactions
lbs <- as.data.frame(cbind(recon3D_BIGG[[12]],recon3D_BIGG[[13]],recon3D_BIGG[[16]]))
lbs$direction <- ifelse((recon3D_BIGG[[13]] + recon3D_BIGG[[12]]) >= 0,"forward","backward")
reversible <- ifelse(recon3D_BIGG[[16]] == 1, TRUE, FALSE)

#get the reaction ids
reaction_ids <- unlist(recon3D_BIGG[[8]])

gene_mapping <- recon2_redhuman$gene_mapping
gene_mapping_vec <- gene_mapping$name
names(gene_mapping_vec) <- as.character(gene_mapping$X1)



#create a dataframe to map reaction indexes with genes
reaction_to_genes <- list()
for(i in 1:length(reaction_list))
{
  if(length(reaction_list[[i]][[1]] > 0))
  {
    lol <- F
    genes <- unique(gsub(" and ","_",gsub("[()]","",gsub("_AT[0-9]+","",strsplit(reaction_list[[i]][[1]], split = " or ")[[1]]))))
    df <- as.data.frame(matrix(NA,length(genes), 2))
    df[,1] <- i
    if(sum(as.character(genes) %in% names(gene_mapping_vec)) > 0)
    {
      for(k in 1:length(as.character(genes)))
      {
        if(as.character(genes[k]) %in% names(gene_mapping_vec))
        {
          genes[k] <- gene_mapping_vec[as.character(genes[k])]
        }
      }
      
    } else
    {
      for(k in 1:length(genes))
      {
        if(grepl("_",genes[k]))
        {
          if(grepl("4967",genes[k]))
          {
            print(genes[k])
            lol <- T
          }
          genes_splitted <- strsplit(genes[k], split = "_")[[1]]
          for(j in 1:length(genes_splitted))
          {
            if(lol) print(genes_splitted[j])
            if(as.character(genes_splitted[j]) %in% names(gene_mapping_vec))
            {
              genes_splitted[j] <- gene_mapping_vec[as.character(genes_splitted[j])]
            }
            if(lol) print(genes_splitted[j])
          }
          genes[k] <- paste0(genes_splitted, collapse = "_")
          if(lol) print(genes)
          if(lol) print(i)
        }
      }
      
    }
    df[,2] <- genes
    reaction_to_genes[[i]] <- df
  } else
  {
    df <- as.data.frame(matrix(NA,1, 2))
    df[,1] <- i
    df[,2] <- reaction_ids[i]
    reaction_to_genes[[i]] <- df
  }
}
reaction_to_genes <- as.data.frame(do.call(rbind,reaction_to_genes))
reaction_to_genes[grepl("[a-z]",reaction_to_genes$V2),"V2"] <- paste("orphanReac",reaction_to_genes[grepl("[a-z]",reaction_to_genes$V2),"V2"],sep = "")
reaction_to_genes <- unique(reaction_to_genes)

reaction_to_genes_original <- reaction_to_genes

#get metabolites
metabolites <- unlist(recon3D_BIGG[[1]])
metabolites_names <- unlist(recon3D_BIGG[[2]])

metabolites_paper <- unlist(recon3D_paper[[15]])

HMDB_paper <- unlist(sapply(recon3D_paper[[16]], function(x){
  if(dim(x[[1]])[1] > 0)
  {
    return(x)
  } else
  {
    return("unknown")
  }
}))

#we get hmdb from the paper version of recon3d model, which has better hmdb mapping than the bigg model one.
HMDB_paper <- as.data.frame(cbind(HMDB_paper, unlist(recon3D_paper[[2]])))
names(HMDB_paper) <- c("HMDB","BIGG")
HMDB_paper$BIGG <- gsub("[[]","_",HMDB_paper$BIGG)
HMDB_paper$BIGG <- gsub("[]]","",HMDB_paper$BIGG)
HMDB_paper[grepl("HMDB[0-9]{5}$",HMDB_paper$HMDB),"HMDB"] <- gsub("HMDB","HMDB00",HMDB_paper[grepl("HMDB[0-9]{5}$",HMDB_paper$HMDB),"HMDB"])
# HMDB_paper$HMDB <- gsub("HMDB","HMDB00",HMDB_paper$HMDB)

##some manual correction
wrong_HMDB_to_right_HMDB <- as.data.frame(
  read_csv("support/wrong_HMDB_to_right_HMDB.csv"))
wrong_HMDB_to_right_HMDB_vec <- wrong_HMDB_to_right_HMDB$right_id
names(wrong_HMDB_to_right_HMDB_vec) <- wrong_HMDB_to_right_HMDB$wrong_id

HMDB_paper[which(HMDB_paper$HMDB %in% names(wrong_HMDB_to_right_HMDB_vec)),"HMDB"] <- wrong_HMDB_to_right_HMDB_vec[HMDB_paper[which(HMDB_paper$HMDB %in% names(wrong_HMDB_to_right_HMDB_vec)),"HMDB"]]

HMDB_paper$HMDB_comp <- ifelse(HMDB_paper$HMDB == "unknown",
                               HMDB_paper$BIGG,
                               paste(HMDB_paper$HMDB, gsub(".*_","",HMDB_paper$BIGG), sep = "_"))

#we complete the paper version mapping with bigg when not available
BIGG_mapping<- as.data.frame(read_delim("support/bigg_iduniversal_bigg_idnamemodel_listdatabase", 
              delim = "\t", escape_double = FALSE, 
              trim_ws = TRUE))

BIGG_mapping <- BIGG_mapping[grepl("HMDB",BIGG_mapping$database_links),]
BIGG_mapping$database_links <- gsub(".*hmdb[/]","",BIGG_mapping$database_links)
BIGG_mapping$database_links <- gsub(";.*","",BIGG_mapping$database_links)
BIGG_mapping$HMDB_comp <- BIGG_mapping$database_links
BIGG_mapping[grepl("HMDB[0-9]{5}$",BIGG_mapping$HMDB_comp),"HMDB_comp"] <- gsub("HMDB","HMDB00",BIGG_mapping[grepl("HMDB[0-9]{5}$",BIGG_mapping$HMDB_comp),"HMDB_comp"])

BIGG_mapping[which(BIGG_mapping$HMDB_comp %in% names(wrong_HMDB_to_right_HMDB_vec)),"HMDB_comp"] <- wrong_HMDB_to_right_HMDB_vec[BIGG_mapping[which(BIGG_mapping$HMDB_comp %in% names(wrong_HMDB_to_right_HMDB_vec)),"HMDB_comp"]]

BIGG_mapping$HMDB_comp <- paste(BIGG_mapping$HMDB_comp,gsub(".*_","",BIGG_mapping$bigg_id),sep = "_")
# BIGG_mapping$HMDB_comp <-gsub("HMDB","HMDB00",BIGG_mapping$HMDB_comp)
BIGG_mapping$bigg_id <- gsub("__","_",BIGG_mapping$bigg_id)

mapping_vec_BIGG <- BIGG_mapping$HMDB_comp
names(mapping_vec_BIGG) <- BIGG_mapping$bigg_id


HMDB_paper$HMDB_comp[which(HMDB_paper$HMDB_comp %in% names(mapping_vec_BIGG))] <- mapping_vec_BIGG[HMDB_paper$HMDB_comp[which(HMDB_paper$HMDB_comp %in% names(mapping_vec_BIGG))]]

mapping_table <- HMDB_paper
mapping_table$name <- unlist(sapply(recon3D_paper[[15]], function(x){
  if(dim(x[[1]])[1] > 0)
  {
    return(x)
  } else
  {
    return("unknown")
  }
}))

mapping_table$formula <- unlist(sapply(recon3D_paper[[13]], function(x){
  if(dim(x[[1]])[1] > 0)
  {
    return(x)
  } else
  {
    return("unknown")
  }
}))

mapping_table$formula_det <- unlist(sapply(recon3D_paper[[14]], function(x){
  if(dim(x[[1]])[1] > 0)
  {
    return(x)
  } else
  {
    return("unknown")
  }
}))

mapping_table$inchi <- unlist(sapply(recon3D_paper[[17]], function(x){
  if(dim(x[[1]])[1] > 0)
  {
    return(x)
  } else
  {
    return("unknown")
  }
}))

save(mapping_table, file = "support/mapping_table.RData")

mapping_vec <- HMDB_paper$HMDB_comp
names(mapping_vec) <- HMDB_paper$BIGG

metabolites <- gsub("__","_",metabolites) # to make it coherent with the names in paper version

metabolites[which(metabolites %in% names(mapping_vec))] <- mapping_vec[metabolites[which(metabolites %in% names(mapping_vec))]]

#We will modify the name of genes as we go so we save the original names
reaction_to_genes <- reaction_to_genes_original

#create the 2 coulmn format network between metabolites and enzymes
reactions_df <- list()

#we do reactions 1 by 1 (colunms of stochiometric matrix)
for(i in 1:length(S[1,]))
{
  print(i)
  #get the reactions stochiometry
  reaction <- S[,i]
  
  #modify gene name so reactions that are catalised by same enzyme stay separated
  reaction_to_genes[reaction_to_genes$V1 == i,2] <- paste(paste("Gene",i, sep = ""),reaction_to_genes[reaction_to_genes$V1 == i,2], sep = "__")
  
  #get the enzymes associated with reaction
  genes <- reaction_to_genes[reaction_to_genes$V1 == i,2]
  
  if(lbs[i,4] == "forward")
  {
    #get reactant
    reactants <- metabolites[reaction == -1]
    
    #get the products
    products <- metabolites[reaction == 1]
  } else
  {
    #get reactant
    reactants <- metabolites[reaction == 1]
    
    #get the products
    products <- metabolites[reaction == -1]
  }
  
  # reactants[which(reactants %in% names(mapping_vec))] <- mapping_vec[reactants[which(reactants %in% names(mapping_vec))]]
  # products[which(products %in% names(mapping_vec))] <- mapping_vec[products[which(products %in% names(mapping_vec))]]
  # 
  reactants <- paste("Metab__",reactants,sep = "")
  products <- paste("Metab__",products,sep = "")
  
  #check how many rows(interactions) will be necessary in the 2 column format of this reaction
  number_of_interations <- length(reactants) + length(products)
  
  #now for each enzymes, we create a two column datframe recapitulating the interactions between the metabolites and this enzyme
  reaction_df <- list()
  j <- 1
  for(gene in genes)
  {
    gene_df <- as.data.frame(matrix(NA,number_of_interations,2))
    gene_df$V1 <- c(reactants,rep(gene,number_of_interations-length(reactants))) #reactants followed by the enzyme (the enzyme is repeated asmany time as they are products)
    gene_df$V2 <- c(rep(gene,number_of_interations-length(products)),products) #enzyme(repeated asmany time as they are reactants) followed by products
    
    if(length(reactants) >= 2 & length(products) >= 2 & length(products) == length(reactants))
    {
      if(sum(gsub("_[a-z]$","",reactants) == gsub("_[a-z]$","",products)) == length(reactants))
      {
        prefix <- str_extract(gene,"Gene[0-9]+__")
        prefix <- gsub("__","",prefix)
        gene <- gsub("Gene[0-9]+__","",gene)
        for(k in 1:length(reactants))
        {
          new_prefix <- paste("0000",k, sep = "")
          new_prefix <- paste(new_prefix,"__",sep = "")
          new_prefix <- paste(prefix,new_prefix,sep = "")
          
          # if(grepl("_e$",gene_df[k,1]) |  grepl("_e$",gene_df[k+length(reactants),2])) #Remove external import and exports
          # {
          #   gene_df[k,] <- c(NA,NA)
          #   gene_df[k+length(reactants),] <- c(NA,NA)
          # } else
          # {
          #   gene_df[k,2] <- paste(new_prefix, gene, sep = "")
          #   gene_df[k+length(reactants),1] <- paste(new_prefix, gene, sep = "")
          # }
          
          gene_df[k,2] <- paste(new_prefix, gene, sep = "")
          gene_df[k+length(reactants),1] <- paste(new_prefix, gene, sep = "")
        }
      }
    }
    
    if(reversible[i]) #if reaction is reversible, we do the same but inverse reactant and products
    {
      gene_df_reverse <- as.data.frame(matrix(NA,number_of_interations,2))
      gene_df_reverse$V1 <- c(rep(paste(gene,"_reverse",sep = ""),number_of_interations-length(products)),products)
      gene_df_reverse$V2 <- c(reactants,rep(paste(gene,"_reverse",sep = ""),number_of_interations-length(reactants)))
      gene_df <- as.data.frame(rbind(gene_df,gene_df_reverse))
    }
    
    # if(sum(grepl("_e$",as.character(c(gene_df[,1],gene_df[,2]))) > 0)) #remove external imports exports
    # {
    #   reaction_df[[j]] <- as.data.frame(matrix(NA,1,2))
    # } else
    # {
    #   reaction_df[[j]] <- gene_df 
    # }
    
    reaction_df[[j]] <- gene_df
    j <- j+1
  }
  #the individual enzyme dataframes of this reaction are combined into one reaction dataframe
  reaction_df <- as.data.frame(do.call(rbind,reaction_df))
  
  #the reaction dataframe is added to the list of all reaction reaction dataframes
  reactions_df[[i]] <- reaction_df
}

#the individual reaction dataframesare combined into one
reactions_df <- as.data.frame(do.call(rbind,reactions_df))

reactions_df <- reactions_df[reactions_df$V1 != "Metab__" & reactions_df$V2 != "Metab__",]

reactions_df <- reactions_df[complete.cases(reactions_df),]

# reactions_df <- reactions_df[!(grepl("Metab__.*_e$",reactions_df$V1)) & !(grepl("Metab__.*_e$",reactions_df$V2)),] # Remove external metabolism
###############################
###############################
###############################

gene_metab_network <- reactions_df
gene_metab_network$V1 <- gsub("_[a-z]$","",gene_metab_network$V1)
gene_metab_network$V2 <- gsub("_[a-z]$","",gene_metab_network$V2)

nodes <- rep(1,length(unique(c(gene_metab_network$V1, gene_metab_network$V2))))
names(nodes) <- unique(c(gene_metab_network$V1, gene_metab_network$V2))

metabs <- nodes[grepl("Metab__",names(nodes))]

i <- 1
for(metab in names(metabs))
{
  print(i)
  metabs[metab] <- sum(metab == gene_metab_network$V1)
  metabs[metab] <- metabs[metab] + sum(metab == gene_metab_network$V2)
  i <- i+1
}

metabs_sorted <- sort(metabs, decreasing = T)

cofactor <- metabs_sorted[metabs_sorted >= 400]
cofactor

metabs_sorted <- metabs_sorted[metabs_sorted < 400] #smallest number of connections before we find important metabolties
#glycine and glutamate both have more than 300 connections lol

reactions_df_no_cofac <- reactions_df[gsub("_[a-z]$","",reactions_df$V1) %in% names(metabs_sorted) | 
                                        gsub("_[a-z]$","",reactions_df$V2) %in% names(metabs_sorted),]

names(reactions_df_no_cofac) <- c("source","target")

reactions_df_no_cofac$source <- gsub("__NA$","__orphanReacNA",reactions_df_no_cofac$source)
reactions_df_no_cofac$target <- gsub("__NA$","__orphanReacNA",reactions_df_no_cofac$target)


write_csv(reactions_df_no_cofac, file = "result/recon3D_bigg_no_cofactor.csv")
