# Dependencies

library(reshape2)
library(rstatix)
library(tidyr)
library(ggpubr)
#library(plyr)
library(dplyr)

source("UDFs/matchchrom_imodes.R")

#load("line260.RData")
#load("Saved/210715_nets.Rdata")

## input and arrange data
neg_data <- read.csv("00_DATA_selected_neg.csv", row.names = 1)
pos_data <- read.csv("00_DATA_selected_pos.csv", row.names = 1)
samplelist <- read.csv("00_200615_Samplelist.csv")

#Remove features with high sd on QC (dispersion index)
dindex_neg <- data.frame(Samples = apply(neg_data[,48:97], 1, sd),
                         QCs = apply(neg_data[,94:98], 1, sd))  
dindex_neg$dindex <- dindex_neg$QCs/dindex_neg$Samples

dindex_pos <- data.frame(Samples = apply(pos_data[,48:97], 1, sd),
                         QCs = apply(pos_data[,94:98], 1, sd))  
dindex_pos$dindex <- dindex_pos$QCs/dindex_pos$Samples

neg_data <- neg_data[rownames(neg_data) %in% rownames(dindex_neg[dindex_neg$dindex < 0.3,]),]
pos_data <- pos_data[rownames(pos_data) %in% rownames(dindex_pos[dindex_pos$dindex < 0.3,]),]
rm(dindex_neg, dindex_pos)

## ID feature ionization modes to combine datasets

rownames(neg_data) <- paste0(rownames(neg_data), "_neg")
rownames(pos_data) <- paste0(rownames(pos_data), "_pos")

neg_feat <- neg_data[,1:47]
neg_feat$ionmode <- "negative"
neg_data <- neg_data[,48:ncol(neg_data)]

pos_feat <- pos_data[,1:47]
pos_feat$ionmode <- "positive"
pos_data <- pos_data[,48:ncol(pos_data)]

colnames(neg_data) <- samplelist$Sample_ID[match(colnames(neg_data), paste0("X", samplelist$filename_LC_neg))]
colnames(pos_data) <- samplelist$Sample_ID[match(colnames(pos_data), paste0("X", samplelist$filename_LC_pos))]

neg_data <- neg_data[,order(colnames(neg_data))]
pos_data <- pos_data[,order(colnames(pos_data))]

## remove QCs - It was already normalized in metaboanalystR
QC_neg <- neg_data[,grep(pattern = "QC", x = colnames(neg_data))]
neg_data[,grep(pattern = "QC", x = colnames(neg_data))] <- NULL

QC_pos <- pos_data[,grep(pattern = "QC", x = colnames(pos_data))]
pos_data[,grep(pattern = "QC", x = colnames(pos_data))] <- NULL

# remove features without MS2 information
neg_feat <- neg_feat[neg_feat$msdial == T,]
neg_data <- neg_data[rownames(neg_data) %in% rownames(neg_feat),]

pos_feat <- pos_feat[pos_feat$msdial == T,]
pos_data <- pos_data[rownames(pos_data) %in% rownames(pos_feat),]

## match neg/pos
neg_feat$Peak_ID <- rownames(neg_feat)
pos_feat$Peak_ID <- rownames(pos_feat)

match_tab <- matchchrom_imodes(peak_table_1 = neg_feat, peak_table_2 = pos_feat,
                          data_table_1 = neg_data, data_table_2 = pos_data,
                          cor_min = 0.8, rt_dev = 0.06, met_annot = T)
colnames(match_tab) <- c(paste0(colnames(match_tab)[1:48], "_neg"),
                         paste0(colnames(match_tab)[49:96], "_pos"))

## Combine data
comb_data <-  merge(t(neg_data), t(pos_data), by = "row.names")
rownames(comb_data) <- comb_data[,1]
comb_data$Row.names <- NULL

## combine feat table
pos_feat[,c(7,8)] <- pos_feat[,c(8, 7)]
colnames(pos_feat)[7:8] <- c("rt_tpj", "mz_tpj")
colnames(neg_feat) == colnames(pos_feat)
comb_feat <- rbind(neg_feat, pos_feat)

# keep table with all
comb_data_all <- as.data.frame(t(comb_data))
comb_feat_all <- comb_feat

## Remove matches btw neg and pos - !Matches are not unique!
pos_data_unique <- pos_data[!rownames(pos_data) %in% match_tab$Peak_ID_pos,]
neg_data_unique <- neg_data[!rownames(neg_data) %in% match_tab$Peak_ID_neg,]

comb_data <- comb_data[,!colnames(comb_data) %in% match_tab$Peak_ID_pos]

comb_feat <- comb_feat[rownames(comb_feat) %in% colnames(comb_data),]
comb_data <- as.data.frame(t(comb_data))

comb_feat$ionmode[comb_feat$Peak_ID %in% match_tab$Peak_ID_neg] <- "Both" 
table(comb_feat$ionmode)
table(comb_feat_all$ionmode)

## average replicates - clustering
avg_comb <- data.frame(row.names = colnames(comb_data),
                      Sample = as.factor(sub("(_.*?)_.*", "\\1", colnames(comb_data))),t(comb_data)
                      )

avg_comb <- avg_comb %>% group_by(Sample) %>% summarise_at(colnames(avg_comb)[-1], mean)
avg_comb <- as.data.frame(avg_comb)
rownames(avg_comb) <- avg_comb$Sample
avg_comb$Sample <- NULL

avg_all <- data.frame(row.names = colnames(comb_data_all),
                      Sample = as.factor(sub("(_.*?)_.*", "\\1", colnames(comb_data_all))),t(comb_data_all)
                      )

avg_all <- avg_all %>% group_by(Sample) %>% summarise_at(colnames(avg_all)[-1], mean)
avg_all <- as.data.frame(avg_all)
rownames(avg_all) <- avg_all$Sample
avg_all$Sample <- NULL

# clustering
dist_mat <- dist(scale(avg_comb), method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)

###### Information theory
temp <- comb_data_all #with all data seperation is better than with match

data_summary <- data.frame(row.names = colnames(temp),
                           Genepool = as.factor(gsub("(.*)_.*", "\\1", colnames(temp))),
                           t((temp))
)

############################### for avgs only

#data_summary <- melt(data_summary, id.vars = "Genepool",
#                     variable.name = "Feature")


#data_summary <- data_summary %>% group_by(Genepool, Feature) %>%
#  summarise(Mean = mean(value, na.rm = TRUE), 
#            Median = median(value, na.rm = TRUE),
#            Sum = sum(value, na.rm = TRUE),
#           SD = sd(value, na.rm = TRUE),
#            CV = sd(value, na.rm = TRUE)/mean(value, na.rm = TRUE),
#            rMAD = mad(value, constant = 1, na.rm = TRUE)/median(value, na.rm = TRUE),
#            NNa = sum(is.na(value)),
#            N = n() - NNa
#  ) %>%
#  as.data.frame()


#Pij_tab <- data_summary %>% pivot_wider(id_cols = c("Feature"),
#                                        names_from = "Genepool", values_from = "Sum") %>%
#  as.data.frame()

#rownames(Pij_tab) <- Pij_tab$Feature
#Pij_tab$Feature <- NULL

###############################

############################### for all run from here
Pij_tab <- as.data.frame(t(data_summary[,-1]))
###############################

Pij_tab <- as.data.frame(t(apply(Pij_tab, 1, function(x) x/colSums(Pij_tab))))

Pij_tab$Pi <- rowMeans(Pij_tab)
Si_tab <- Pij_tab/Pij_tab$Pi
Si_tab$Pi <- NULL
Si_tab <- Si_tab*log2(Si_tab)
Si_tab$Si <- rowSums(Si_tab)/ncol(Si_tab)

tj <- Pij_tab*Si_tab$Si
H_tab <- Pij_tab*log2(Pij_tab)

Infothe <- data.frame(row.names = names(H_tab),
                      Hj = -colSums(H_tab, na.rm = T),
                      tj = colSums(tj, na.rm = T))
Infothe <- Infothe[-nrow(Infothe),]
temp <- as.factor(sub("_.*", "", rownames(Infothe)))

Infothe <- rbind(Infothe,
      t(data.frame(MW = colMeans(Infothe[temp == "MW",-3]),
                   MD = colMeans(Infothe[temp == "MD",-3]))))

Infothe$Genepool <- as.factor(sub("_.*", "", rownames(Infothe)))


p <- ggplot(Infothe, aes(x=Hj, y=tj, label=rownames(Infothe), color = Genepool)) + geom_text() + 
  xlab("Hj (Diversity)") + ylab("\u03B4j (Specialization)") + theme_minimal() +
  scale_color_manual(labels = c("MW", "MD"), values = c("#d8b365", "#5ab4ac"))
p


##################### for avg only
#md <- Infothe[1:7,]
#mw <- Infothe[8:15,]
#####################


##################### for all
md <-  Infothe[1:20,]
mw <- Infothe[21:46,]
#####################

# Test
var.test(md$Hj, mw$Hj)
shapiro.test(md$Hj)
shapiro.test(mw$Hj)
t.test(md$Hj, mw$Hj, var.equal = T)

var.test(md$tj, mw$tj)
shapiro.test(md$tj)
shapiro.test(mw$tj)
t.test(md$tj, mw$tj, var.equal = T)

rm(temp)

##########################################################################


###### Univariate stats

# Summary statistics per genepool

source("UDFs/Summ_data_bean_2.R")
source("UDFs/stat_bean.R")

data_summary <- Summ_data_bean(log(comb_data))
data_summary_all <- Summ_data_bean(log(comb_data_all))

data_summary <- stat_bean(log(comb_data), data_summary)
data_summary_all <- stat_bean(log(comb_data_all), data_summary_all)

#write.csv(data_summary, "00_DATA_SummaryTable.csv")

# Add info theory information to table 
temp <- data.frame(row.names = rownames(Si_tab),
           Si = Si_tab$Si)
data_summary <- merge(data_summary, temp, by = "row.names")
rownames(data_summary) <- data_summary$Row.names
data_summary$Row.names <- NULL

data_summary_all <- merge(data_summary_all, temp, by = "row.names")
rownames(data_summary_all) <- data_summary_all$Row.names
data_summary_all$Row.names <- NULL

# Add proportions of in MW and MD to table  
temp <- exp(data_summary[,c(2,8)])/rowSums(exp(data_summary[,c(2,8)]))
colnames(temp) <- c("Prop_MW", "Prop_MD")
data_summary <- cbind(data_summary, temp)

temp <- exp(data_summary_all[,c(2,8)])/rowSums(exp(data_summary_all[,c(2,8)]))
colnames(temp) <- c("Prop_MW", "Prop_MD")
data_summary_all <- cbind(data_summary_all, temp)


##########################################################################

########################## Check Distribution of mean
# All data
mean_plot <- data_summary[, c(1, 7)] 

mean_plot <- melt(mean_plot)

mean_plot %>% wilcox_test(value ~ variable)
mu <- mean_plot %>% group_by(variable) %>% summarise(grp.mean=mean(value)) %>% as.data.frame()

p <- ggplot(mean_plot, aes(value, color = variable)) + geom_density() +
  theme_minimal() + labs(x = "mean", y="Density", color = "Genepool") +
  scale_color_manual(labels = c("MW", "MD"), values = c("#d8b365", "#5ab4ac")) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),
             linetype="dashed")
p
#ggsave(filename = "mean_plot.pdf", useDingbats=FALSE,
#       plot = p)


########################## Check Distribution of CV
# only significant Bonferoni
#CV_plot <- data_summary[data_summary$p.adj.signif_bonferoni != "ns", c(4, 10)] 

# All data
CV_plot <- data_summary[ c(4, 10)] 

CV_plot <- melt(CV_plot)

CV_plot %>% wilcox_test(value ~ variable)
mu <- CV_plot %>% group_by(variable) %>% summarise(grp.mean=mean(value)) %>% as.data.frame()

p <- ggplot(CV_plot, aes(value, color = variable)) + geom_density() +
  theme_minimal() + labs(x = "CV", y="Density", color = "Genepool") +
  scale_color_manual(labels = c("MW", "MD"), values = c("#d8b365", "#5ab4ac")) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),
             linetype="dashed")

p
#ggsave(filename = "Figure 2D.pdf", useDingbats=FALSE,
#       plot = p)

# CV analysis

CV_summary <- data.frame(row.names = c("DEM_FDR", "NDEM_FDR", "DEM_BON", "NDEM_BON", "TOTAL"),
                         CV_MW = c(
                           mean(data_summary$CV_MW[data_summary$p.adj.signif_fdr != "ns"]),
                           mean(data_summary$CV_MW[data_summary$p.adj.signif_fdr == "ns"]),
                           mean(data_summary$CV_MW[data_summary$p.adj.signif_bonferoni != "ns"]),
                           mean(data_summary$CV_MW[data_summary$p.adj.signif_bonferoni == "ns"]),
                           mean(data_summary$CV_MW)
                         ),
                         CV_MD = c(
                           mean(data_summary$CV_MD[data_summary$p.adj.signif_fdr != "ns"]),
                           mean(data_summary$CV_MD[data_summary$p.adj.signif_fdr == "ns"]),
                           mean(data_summary$CV_MD[data_summary$p.adj.signif_bonferoni != "ns"]),
                           mean(data_summary$CV_MD[data_summary$p.adj.signif_bonferoni == "ns"]),
                           mean(data_summary$CV_MD)
                         ))

CV_summary$L_CV <- 1-(CV_summary$CV_MD/CV_summary$CV_MW)
#write.csv(CV_summary, "CV_summary.csv")

##########################

# Check distribution of metabolite specificty 
# according to genepools predominance

Si_dist <- data.frame(variable = "MW", value = data_summary$Si[data_summary$Prop_MW > 0.5])
temp <- data.frame(variable = "MD", value = data_summary$Si[data_summary$Prop_MW < 0.5])

Si_dist <- rbind(Si_dist, temp)

Si_dist %>% wilcox_test(value ~ variable)
mu <- Si_dist %>% group_by(variable) %>% summarise(grp.mean=mean(value)) %>% as.data.frame()

p <- ggplot(Si_dist, aes(value, color = variable)) + geom_density() +
  theme_minimal() + labs(x = "Si", y="Density", color = "Genepool") +
  scale_color_manual(labels = c("MD", "MW"), values = c( "#5ab4ac", "#d8b365")) +
  geom_vline(data=mu, aes(xintercept=grp.mean, color=variable),
             linetype="dashed")

p
#ggsave(filename = "Figure_Sidensity.pdf", useDingbats=FALSE,
#       plot = p)


###### Networks assembly

### Similarity network
library(OrgMassSpecR)
source("UDFs/spec_vec_to_lib.R")
source("UDFs/flattenCorrMatrix.R")

#convert all MSDIAL MS2 intensity vectors into a list of 2 column dataframes to input into Spec similarity fun
neg_ms2 <- spec_vec_to_lib(neg_feat$MS.MS.spectrum_msdial)

#Apply spec sim calculation to all pairs in the list 
simnet_neg <- sapply(neg_ms2, function(x) sapply(neg_ms2,
                                                  function(y) SpectrumSimilarity(x,y,t = 0.01, print.graphic = F))) # ~ 1 sec per 200 pairs - 1 hour for ~880 features

rm(neg_ms2)

backup_simnet_neg <- simnet_neg


#fix colnames
colnames(simnet_neg) <- neg_feat$Peak_ID
rownames(simnet_neg) <- neg_feat$Peak_ID

# change to long format removing lower triangle and diagonal and Nas
simnet_neg <- flattenCorrMatrix(simnet_neg)

#simnet_neg <- simnet_neg[!is.na(simnet_neg$Similarity_score),]

# filter by similarity score
#simnet_neg <- simnet_neg[simnet_neg$Similarity_score > 0.7,]

#simnet_neg$edge_sim <- "s"

# write csv file for cytoscape
#write.csv(simnet_neg, "out_netsim2.csv")


#convert all MSDIAL MS2 intensity vectors into a list of 2 column dataframes to input into Spec similarity fun
pos_ms2 <- spec_vec_to_lib(pos_feat$MS.MS.spectrum_msdial)

#Apply spec sim calculation to all pairs in the list 
simnet_pos <- sapply(pos_ms2, function(x) sapply(pos_ms2,
                                                 function(y) SpectrumSimilarity(x,y,t = 0.01, print.graphic = F))) # ~ 1 sec per 200 pairs - 1 hour for ~880 features

rm(pos_ms2)

backup_simnet_pos <- simnet_pos


#fix colnames
colnames(simnet_pos) <- pos_feat$Peak_ID
rownames(simnet_pos) <- pos_feat$Peak_ID

# change to long format removing lower triangle and diagonal and Nas
simnet_pos <- flattenCorrMatrix(simnet_pos)
#simnet_pos <- simnet_pos[!is.na(simnet_pos$Similarity_score),]

# filter by similarity score
#simnet_pos <- simnet_pos[simnet_pos$Similarity_score > 0.7,]

#simnet_pos$edge_sim <- "s"



### Correlation network

# create correlation matrix using spearman rank correlation
cornet_neg <- cor(t(neg_data), method = "spearman")

# uniform feature name for all networks
colnames(cornet_neg) <- gsub("\\_.*", "", colnames(cornet_neg))
rownames(cornet_neg) <- gsub("\\_.*", "", rownames(cornet_neg))

# change to long format removing lower triangle and diagonal
cornet_neg <- flattenCorrMatrix(cornet_neg)
colnames(cornet_neg)[3] <- "spearman_cor"

# Filter by correlation coef
#cornet_neg <- cornet_neg[cornet_neg$spearman_cor > 0.7,]

#cornet_neg$edge_cor <- "c"

# write csv file for cytoscape
#write.csv(cornet_neg, "out_netcor.csv")


# create correlation matrix using spearman rank correlation
cornet_pos <- cor(t(pos_data), method = "spearman")

# uniform feature name for all networks
colnames(cornet_pos) <- gsub("\\_.*", "", colnames(cornet_pos))
rownames(cornet_pos) <- gsub("\\_.*", "", rownames(cornet_pos))

# change to long format removing lower triangle and diagonal
cornet_pos <- flattenCorrMatrix(cornet_pos)
colnames(cornet_pos)[3] <- "spearman_cor"

# Filter by correlation coef
#cornet_pos <- cornet_pos[cornet_pos$spearman_cor > 0.7,]

#cornet_pos$edge_cor <- "c"

# write csv file for cytoscape
#write.csv(cornet_pos, "out_netcor_pos.csv")



### Bioreact network
source("UDFs/Net_bio.R")

bionet_neg <- net.bio(neg_feat, bio_loss = read.csv("00_Neutral_losses.csv"), mz_dev = 0.007) # With 0.007 there were no duplicates

#bionet_neg$edge_bio <- "b"

#write.csv(bionet_neg, "out_netbio.csv")


bionet_pos <- net.bio(pos_feat, bio_loss = read.csv("00_Neutral_losses.csv"), mz_dev = 0.007) # With 0.007 there were no duplicates

#bionet_pos$edge_bio <- "b"

#write.csv(bionet_pos, "out_netbio_pos.csv")

# Merge all three networks and export to cytoscape
allnet <- merge(bionet_neg, simnet_neg, by = "row.names", all = T)
row.names(allnet) <- allnet$Row.names
allnet$Row.names <- NULL
allnet <- merge(allnet, cornet_neg, by = "row.names", all = T)
allnet <- allnet[!is.na(allnet$Similarity_score),]

allnet <- allnet[allnet$Similarity_score > 0.7,]
#back_allnetneg <- allnet
allnet$ionmode <- "negative"


temp <- merge(bionet_pos, simnet_pos, by = "row.names", all = T)
row.names(temp) <- temp$Row.names
temp$Row.names <- NULL
temp <- merge(temp, cornet_pos, by = "row.names", all = T)
temp <- temp[!is.na(temp$Similarity_score),]

temp <- temp[temp$Similarity_score > 0.7,]
#back_allnetpos <- temp
temp$ionmode <- "positive"


# merge both ionmodes
allnet <- rbind(allnet, temp)
allnet <- allnet[,-c(1,2,3,12,13)]
#write.csv(allnet, "net_all.csv", row.names = F)

allnet_cor <- allnet[abs(allnet$spearman_cor) > 0.5 &
                       allnet$Feature_A.y %in% rownames(comb_data) &
                       allnet$Feature_B.y %in% rownames(comb_data), ]
#write.csv(allnet_cor, "net_all_cor.csv", row.names = F)


#save.image("210715_nets.Rdata")

# export feature info - have to edit tpj class for some reason otherwise cytoscape can#t read
#write.csv(neg_feat, "out_feat_neg.csv", row.names = F)



##############################################################################
# Get 1st PC for each module
#modules_data <- read.csv("00_DATA_net_modules_bigger2.csv", row.names = 1)
#modules_data[,-1] <- log(modules_data[,-1])

modules_data <- read.csv("00_DATA_Modules_mergedNets.csv", row.names = 1)
modules_data <- merge(modules_data, t(avg_all), by = "row.names")
rownames(modules_data) <- modules_data$Row.names
modules_data$Row.names <- NULL
modules_data[,-1] <- log(modules_data[,-1])

for(i in 1:max(modules_data$Module)){
  if(i == 1){
    
    module_group <- modules_data[modules_data$Module == i,-1]
    pca_module <- prcomp(t(module_group), center = T, scale. = T)
    data_PCs <- data.frame(Module_1 = pca_module$x[,1])
  } else{
    
    module_group <- modules_data[modules_data$Module == i,-1]
    pca_module <- prcomp(t(module_group), center = T, scale. = T)
    data_PCs <- cbind(data_PCs, data.frame(pca_module$x[,1]))
    colnames(data_PCs)[ncol(data_PCs)] <- paste0("Module_", i)
  }
}

library(ggbiplot)
pca_comb <- prcomp(data_PCs, center = T, scale. = T)
ggbiplot(pca_comb, var.axes=FALSE, ellipse = T, groups=gsub("_.*", "", rownames(data_PCs)),
         labels = rownames(data_PCs)) +
  ggtitle("") + theme_minimal()


# clustering
dist_mat <- dist(data_PCs, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)


library(pheatmap)
M <- cor(data_PCs)
pheatmap(M, show_rownames = T, show_colnames = T, fontsize = 6)

#write.csv(t(data_PCs), "out_DATA_ModulePC_merged.csv")

### Export files to Sirius

for(i in 1:max(data_summary_all$Modules, na.rm = T)){
  feature_data <- comb_feat_all[rownames(comb_feat_all) %in% 
                                  rownames(data_summary_all[data_summary_all$Modules == i,]),]
  
  source("UDFs/tab_to_sirius.R")
  tab_to_sirius(feature_data, filename = paste0("Sirius_Module_", i))
}


##################################################################

##################################################################
#run last, ggbiplot conflict with dplyr

# PCA of filtered data
theme_set(theme_minimal() + theme(legend.position = "bottom"))

pca_comb <- prcomp(log(avg_comb), center = T, scale. = T)

p <- ggbiplot::ggbiplot(pca_comb, var.axes=FALSE, ellipse = T, groups=gsub("_.*", "", rownames(avg_comb)),
                        labels = rownames(avg_comb)) +
  ggtitle("PCA of LCMS neg+pos ions - log data avgs")
p


pca_comb <- prcomp(log(t(comb_data)), center = T, scale. = T)

p <- ggbiplot::ggbiplot(pca_comb, var.axes=FALSE, ellipse = T, groups=gsub("_.*", "", rownames(t(comb_data))),
                        labels = rownames(t(comb_data))) +
  ggtitle("PCA of LCMS neg+pos ions - log data")
p

pca_varcontribution <- as.data.frame(pca_comb$rotation[,1:2])
head(pca_varcontribution[order(abs(pca_varcontribution$PC1), decreasing = T),])

write.csv(pca_varcontribution, "01_pca_varcontribution.csv")
##################################################################

##########################################################################
# summarize modules
# input data from Networks!
modules <- read.csv("00_DATA_Modules_mergedNets.csv", row.names = 1)
data_summary <- merge(data_summary, modules, by = "row.names", all.x = T)
rownames(data_summary) <- data_summary$Row.names
data_summary$Row.names <- NULL

data_summary_all <- merge(data_summary_all, modules, by = "row.names", all.x = T)
rownames(data_summary_all) <- data_summary_all$Row.names
data_summary_all$Row.names <- NULL

module_summary <- data_summary %>% group_by(Modules) %>%
  summarise(Mean_MWprop = mean(Prop_MW, na.rm = TRUE), 
            Median_MWprop = median(Prop_MW, na.rm = TRUE),
            Mean_Si = mean(Si, na.rm = TRUE),
            Meadian_Si = median(Si, na.rm = TRUE),
            n = n()) %>%
  as.data.frame()

#write.csv(module_summary, "Summary_modules.csv")

# 
sel_m <- 12
data_summary$dummy <- "All"
data_summary$dummy[which(data_summary$Modules == sel_m)] <- paste0("M", sel_m)

# Change density plot line colors by groups
ggplot(data_summary, aes(x=Si, color=dummy)) +
  geom_density() + theme_minimal() +
  scale_color_manual(labels = c("All", paste0("M", sel_m)), values = c("#d8b365", "#5ab4ac"))


##########################################################################






