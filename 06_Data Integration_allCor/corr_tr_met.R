

## Read data transcripts
tr_data_raw <- read.csv("00_DATA_RNAseq.csv")

## Select transcripts under selection according to Elena's paper - restricting only to DE PCA looks much betters
tr_data_raw <- tr_data_raw[which(tr_data_raw$DE.NDE..We. == "DE"|tr_data_raw$DE.NDE.no.We.=="DE"),]

rownames(tr_data_raw) <- paste("tr", 1:nrow(tr_data_raw), sep = "_")#### Get gene ids blasting

#### Get gene ids blasting
# write fasta file to blast
#library(seqinr)
#write.fasta(sequences = as.list(tr_data_raw$Sequence),
#            names = rownames(tr_data_raw),
#            file.out = "PV_fasta.csv")

# read in blast results !adjusted with sed in Linux from blast_result.txt to blast_clean.txt
gene_ids <- read.delim("blast_clean.txt", header = F)
gene_ids <- gene_ids[!grepl(">", gene_ids$V1),]

gene_ids <- data.frame(transcript = gene_ids[grep("Query", gene_ids)],
                       gene_id = gene_ids[grep("Query", gene_ids)+1])

gene_ids$gene_id[grep("Query", gene_ids$gene_id)] <- NA

gene_ids$gene_id <- gsub(".*\\(","",gene_ids$gene_id)
gene_ids$gene_id <- gsub("\\).*","",gene_ids$gene_id)

tr_data_raw$gene_id <- gene_ids$gene_id
#write.csv(tr_data_raw, "DE_transcripts.csv")
#### #### #### #### #### 

# fix names
rownames(tr_data_raw) <- paste("tr", 1:nrow(tr_data_raw), sep = "_")
tr_data <- tr_data_raw[,2:18]

# read metabolite modules data
neg_data <- read.csv("out_DATA_ModulePC_merged.csv", row.names = 1)

# removing missing metab samples from transcript data
tr_data <- tr_data[, colnames(tr_data) %in% colnames(neg_data)]

# remove constant variables
tr_data <- tr_data[apply(tr_data, 1, var, na.rm=TRUE) != 0,]

# match order in both datasets
neg_data <- neg_data[, colnames(tr_data)]

# transpose datasets
neg_data <- as.data.frame(t(neg_data))
tr_data <- as.data.frame(t(tr_data))

# correlate data
cor_data <- cor(scale(tr_data), neg_data, method = "spearman")

# change data to long format
cor_data <- reshape2::melt(data = cor_data, value.name = "cor")
names(cor_data)[1:2] <- c("transcript", "met_module")

# select top 50 per module
library(dplyr)
cor_data <- cor_data %>% 
  arrange(desc(abs(cor))) %>% 
  group_by(met_module) %>% slice(1:50)


# add relevant info
cor_data$gene_id <- tr_data_raw$gene_id[match(cor_data$transcript, rownames(tr_data_raw))]
cor_data$sequence <- tr_data_raw$Sequence[match(cor_data$transcript, rownames(tr_data_raw))]
cor_data$region_soy <- tr_data_raw$blast.con.SOYA[match(cor_data$transcript, rownames(tr_data_raw))]
cor_data$blast_soy <- tr_data_raw$blast.con.SOYA.4[match(cor_data$transcript, rownames(tr_data_raw))]

# add phytozome link
phy <- read.csv("Phytozome_results.csv")
phy <- phy %>% 
  arrange(E.value) %>% 
  group_by(Query.ID) %>% slice(1)

cor_data$phytozome <- phy$View[match(cor_data$transcript, phy$Query.ID)]

# read. manually curated data
add_info <- read.csv("tr_unique_selected_manual.csv")

# add manual info
cor_data$gene_id[which(is.na(cor_data$gene_id))] <- add_info$X[
  match(cor_data$transcript[which(is.na(cor_data$gene_id))], add_info$transcript)]

# fix names from NCBI Blast
cor_data$gene_id <- gsub("PHAVU_", "Phvul.", cor_data$gene_id)
cor_data$gene_id <- gsub("g", "", cor_data$gene_id)

# write file
write.csv(cor_data, "out_corPmodule.csv", row.names = F)



################################
##### Repeat for all trancripts 
################################
## Read data transcripts
tr_data_raw <- read.csv("00_DATA_RNAseq.csv")

# fix names
rownames(tr_data_raw) <- paste("tr", 1:nrow(tr_data_raw), sep = "_")
tr_data <- tr_data_raw[,2:18]

# read metabolite modules data
neg_data <- read.csv("out_DATA_ModulePC_merged.csv", row.names = 1)

# removing missing metab samples from transcript data
tr_data <- tr_data[, colnames(tr_data) %in% colnames(neg_data)]

# remove constant variables
tr_data <- tr_data[apply(tr_data, 1, var, na.rm=TRUE) != 0,]

# match order in both datasets
neg_data <- neg_data[, colnames(tr_data)]

# transpose datasets
neg_data <- as.data.frame(t(neg_data))
tr_data <- as.data.frame(t(tr_data))

# correlate data
cor_data <- cor(scale(tr_data), neg_data, method = "spearman")

# change data to long format
cor_data <- reshape2::melt(data = cor_data, value.name = "cor")
names(cor_data)[1:2] <- c("transcript", "met_module")

# select top 50 per module
library(dplyr)
cor_data <- cor_data %>% 
  arrange(desc(cor)) %>% 
  group_by(met_module) %>% slice(1:50) %>%
  as.data.frame()

cor_data <- cor_data[abs(cor_data$cor) > 0.7,]

sel_tr <- unique(cor_data$transcript)


############################################
# read in blast results !adjusted with sed in Linux from blast_result.txt to blast_clean.txt
gene_ids <- read.delim("Pv_pcor_clean.txt", header = F)

gene_ids <- gene_ids[!grepl(">", gene_ids$V1),]

gene_ids <- data.frame(transcript = gene_ids[grep("Query", gene_ids)],
                       gene_id = gene_ids[grep("Query", gene_ids)+1])

gene_ids$gene_id[grep("Query", gene_ids$gene_id)] <- NA

gene_ids$gene_id <- gsub(".*\\(","",gene_ids$gene_id)
gene_ids$gene_id <- gsub("\\).*","",gene_ids$gene_id)

gene_ids$transcript <- sel_tr[order(sel_tr)]

# fix names from NCBI Blast
gene_ids$gene_id <- gsub("PHAVU_", "Phvul.", gene_ids$gene_id)
gene_ids$gene_id <- gsub("g", "", gene_ids$gene_id)

cor_data$gene_id <- gene_ids$gene_id[match(cor_data$transcript, gene_ids$transcript)]
cor_data$region_soy <- tr_data_raw$blast.con.SOYA[match(cor_data$transcript, rownames(tr_data_raw))]
cor_data$blast_soy <- tr_data_raw$blast.con.SOYA.4[match(cor_data$transcript, rownames(tr_data_raw))]
cor_data$blast_ara <- tr_data_raw$blast.con.Arabidopsis.6[match(cor_data$transcript, rownames(tr_data_raw))]

cor_data <- cor_data[!cor_data$gene_id %in% NA,]

#write.csv(cor_data, "out_DATA_positiveCorModules.csv")









