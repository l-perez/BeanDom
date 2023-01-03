
library(topGO)
library(biomaRt)

#########################################################################
# create list with all genes and GO annotation from Phytozome Pv genome #
#########################################################################

# get access to Phytozome data to get GO terms for gene IDs
phytozome_v13 <- useMart(biomart = "phytozome_mart", 
                         dataset = "phytozome", 
                         host = "https://phytozome-next.jgi.doe.gov")

# get GO terms with gene ids
go_pv <- getBM(attributes = c("gene_name1", "go_id"), 
               filters = "organism_id", 
               values = 442, 
               mart = phytozome_v13)  



# create gene list for topGO
gene2go <- split(go_pv, go_pv$gene_name1)
gene2go <- lapply(gene2go, function(x) unlist(x[2]))

# import selected transcripts from PLS regression
tr_PLS <- read.csv("sel_tr_PLS.csv")
tr_info <- read.csv("tr_unique_selected_manual.csv")

tr_PLS <- tr_info$gene_id[tr_info$transcript %in% tr_PLS$X.3]

# read in selected transcripts per modules
#top_tr <- read.csv("out_corPmodule.csv")
top_tr <- read.csv("out_DATA_positiveCorModules.csv")

######################################################
# Select PLS transcripts and create topGO data 
######################################################

geneNames <- names(gene2go)
geneList <- factor(as.integer(geneNames %in% tr_PLS))
names(geneList) <- geneNames

# run GO analysis ! run with both Bp and MF
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              nodeSize = 5,
              annot = annFUN.gene2GO, gene2GO = gene2go)

######################################################
# Run GO enrichment analysis
######################################################

# run Fisher's test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# result table
tableFisher <- GenTable(GOdata,
                   classicFisher = resultFisher, topNodes = 15)

tableFisher$set <- "PLS"

# save graph with GO terms relationship
printGraph(GOdata, resultFisher, firstSigNodes = 15, fn.prefix = "PLS_tGO", useInfo = "all", pdfSW = TRUE)

#########################################################
### run topGO for all clusters and export graph 
#########################################################

for(n in unique(top_tr$met_module)){
  # create list of selected genes
  geneNames <- names(gene2go)
  geneList <- factor(as.integer(geneNames %in% top_tr$gene_id[top_tr$met_module %in% n]))
  names(geneList) <- geneNames
  
  # run GO analysis
  GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                nodeSize = 5,
                annot = annFUN.gene2GO, gene2GO = gene2go)
  
  # results
  resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  
  ### subgraph induced by the 5 most significant GO terms 
  printGraph(GOdata, resultFisher, firstSigNodes = 15, fn.prefix = n, useInfo = "all", pdfSW = TRUE)
  
  # Append results to table generated for PLS 
  temp <- GenTable(GOdata,
                   classicFisher = resultFisher, topNodes = 15)
  temp$set <- n
  tableFisher <- rbind(tableFisher, temp)
  
  
}

write.csv(tableFisher, "GOresult_all_BP.csv")

#########################################################
### check  
#########################################################

# load results from loop over all modules
tableFisher <- read.csv("GOresult_all.csv")

# Check GO frequency in top significant GOs accross all modules
GO_frequency <- tableFisher[,c(1:4)]
temp <- table(GO_frequency$GO.ID)
GO_frequency$n <- temp[match(GO_frequency$GO.ID, names(temp))]
GO_frequency <- unique(GO_frequency[,2:5])
GO_frequency <- GO_frequency[order(GO_frequency$n, decreasing = T),]

rm(temp)

# crosscheck genes with significant GO
genes <- go_pv$gene_name1[go_pv$go_id %in% GO_frequency$GO.ID[4] & go_pv$gene_name1 %in% tr_info$gene_id]
genes <- go_pv$gene_name1[go_pv$go_id %in% GO_frequency$GO.ID[4]]

temp <- tr_info[tr_info$gene_id %in% genes,]
write.csv(temp, "temp.csv")

