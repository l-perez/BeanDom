# Working directory --------------------------------------------------------------


# Dependencies -------------------------------------------------------------------
library(MetaboAnalystR)
library(dplyr)
library(ggplot2)
theme_set(theme_bw())
library(reshape2)

# Input files -------------------------------------------------------------------
samplelist <- "../00_200615_Samplelist.csv"
datafile <- "00_DATA_bean_metaboanalyst_filterdNN.csv"

# Remove samples/classes (e.g. blanks and technical controls)--------------------------------------------------------
rm_sample <- c("55_neg.mzXML.filtered","57_neg.mzXML.filtered", "60_neg.mzXML.filtered",
               "62_neg.mzXML.filtered","64_neg.mzXML.filtered","5_neg.mzXML.filtered","53_neg.mzXML.filtered")
# "5_neg.mzXMLfiltered","53_neg.mzXMLfiltered"

rm_class <- c("AD")

# end user input ----------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# Timestamp ---------------------------------------------------------------------
startTime <- Sys.time()
startTimeStamp <- format(startTime, "%Y%m%dh%H%M")

# Data input
samplelist = read.csv(samplelist)

mSet<-InitDataObjects("pktable", "stat", FALSE)
mSet<-Read.TextData(mSet, datafile, "colu", "disc")
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);

# Filter variables and samples
#mSet<-FilterVariable(mSet, "sd", "T", 25)
mSet<-PreparePrenormData(mSet)
feature.nm.vec <- c("")
smpl.nm.vec <- rm_sample
grp.nm.vec <- rm_class
mSet<-UpdateData(mSet)

# RLD plot of edited data
p <- RLDplot(samplelist, set="edit", fix_name = ".mzXML.filtered")
ggsave(filename = paste(startTimeStamp, "RLD_rawEdited.svg", sep = "_"))
ggsave(filename = paste(startTimeStamp, "RLD_rawEdited.pdf", sep = "_"))

n <- "GroupPQN"
# Normalization
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, n, "LogNorm", "AutoNorm", ref = "QC", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, paste(startTimeStamp, n, "vnorm_", sep = "_"), "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, paste(startTimeStamp, n, "snorm_", sep = "_"), "png", 72, width=NA)

# Analysis
mSet<-PCA.Anal(mSet)
mSet<-PlotPCA2DScore(mSet, paste(startTimeStamp, n, "pca_score2d_", sep = "_"), "png", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotPCA2DScore(mSet, paste(startTimeStamp, n, "pca_score2d_", sep = "_"), "svg", 72, width=NA, 1,2,0.95,0,0)
mSet<-PlotPCA2DScore(mSet, paste(startTimeStamp, n, "pca_score2d_", sep = "_"), "pdf", 72, width=NA, 1,2,0.95,0,0)

# Normalization without scaling
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, n, "none", "none", ref = "QC", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, paste(startTimeStamp, n, "vnorm_", sep = "_"), "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, paste(startTimeStamp, n, "snorm_", sep = "_"), "png", 72, width=NA)

# prepare output
data_exp <- as.data.frame(cbind(Sample_ID = gsub(pattern = ".mzXML.filtered", "", row.names(mSet[[c("dataSet", "norm")]])),
                            Class = as.character(mSet[[c("dataSet", "edit.cls")]]),
                            as.data.frame(mSet[[c("dataSet", "norm")]])))
write.csv(data_exp, "00_DATA_normalized_PQN.csv")




