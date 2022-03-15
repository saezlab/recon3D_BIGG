
library(readr)
library(reshape2)

HMDB_full_mapping <- as.data.frame(read_csv("support/HMDB_full_mapping.csv")) #from parse_hmdb.ipynb
HMDB_full_mapping <- HMDB_full_mapping[!is.na(HMDB_full_mapping[,2]),]
names(HMDB_full_mapping)[1] <- "good_id"
HMDB_full_mapping <- melt(HMDB_full_mapping, id.vars = "good_id")
HMDB_full_mapping <- HMDB_full_mapping[,c(1,3)]
HMDB_full_mapping <- HMDB_full_mapping[complete.cases(HMDB_full_mapping),]

names(HMDB_full_mapping) <- c("good_HMDB_id","legacy_HMDB_id")

write_csv(HMDB_full_mapping, file = "support/HMDB_full_mapping_long.csv")


