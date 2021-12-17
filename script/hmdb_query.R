library(hmdbQuery)
library(readr)

recon3D_bigg_no_cofactor <- as.data.frame(
  read_csv("result/recon3D_bigg_no_cofactor.csv"))

metabs <- unique(c(recon3D_bigg_no_cofactor$source, recon3D_bigg_no_cofactor$target))
metabs <- metabs[grepl("Metab__",metabs)]
metabs <- metabs[grepl("HMDB",metabs)]
metabs <- gsub("Metab__","",metabs)
metabs <- unique(gsub("_[a-z]$","",metabs))

metabs_df <- as.data.frame(metabs)
metabs_df$super_class <- NA
metabs_df$class <- NA
metabs_df$sub_class <- NA

for(i in 1:length(metabs))
{
  print(paste(i,metabs[i],sep = " "))
  hmdb_query <- tryCatch(
    {

      message("This is the 'try' part")
      hmdbQuery::HmdbEntry(id = metabs[i])
      
    },
    error=function(cond) {
      message(cond)
      
      return(hmdb_query)
    }
  ) 
  
  hmdb_query <- tryCatch(
    {
      
      message("This is the 'try' part")
      hmdbQuery::store(hmdb_query)
      
    },
    error=function(cond) {
      message(cond)
      
      return(hmdb_query)
    }
  )  
  out <- tryCatch(
    {
      
      message("This is the 'try' part")
      metabs_df[i,"super_class"] <- hmdb_query$taxonomy$super_class
      metabs_df[i,"class"] <- hmdb_query$taxonomy$class
      metabs_df[i,"sub_class"] <- hmdb_query$taxonomy$sub_class
      
    },
    error=function(cond) {
      message(cond)
      
      return(NA)
    }
  )  
}

save(metabs_df, file = "Dropbox/recon3D_BIGG/support/HMDB_metab_class.RData")
