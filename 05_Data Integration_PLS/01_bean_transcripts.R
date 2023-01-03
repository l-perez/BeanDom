


## Read data
tr_data_raw <- read.csv("00_DATA_RNAseq.csv")
## Select transcripts under selection according to plant cell paper
tr_data_raw <- tr_data_raw[which(tr_data_raw$DE.NDE..We. == "DE"|tr_data_raw$DE.NDE.no.We.=="DE"),]

rownames(tr_data_raw) <- paste("tr", 1:nrow(tr_data_raw), sep = "_")
tr_data <- tr_data_raw[,2:18]


# read metabolites data
#neg_data <- as.data.frame(t(read.csv("00_DATA_avg_comb.csv", row.names = 1)))
neg_data <- read.csv("out_DATA_ModulePC_bigger2.csv", row.names = 1)

# removing missing metab samples from transcript data
tr_data <- tr_data[, colnames(tr_data) %in% colnames(neg_data)]

# remove constant variables
tr_data <- tr_data[apply(tr_data, 1, var, na.rm=TRUE) != 0,]

# order
neg_data <- neg_data[, colnames(tr_data)]

# sample classes
samples <- gsub("_.*", "", colnames(neg_data))

########### PCA transcripts
# prepare data for PCA
pca_data <- t(tr_data)
pca_data[is.na(pca_data)] <- min(pca_data, na.rm = T)
pca_data[is.na(pca_data)]

# remove constant variables
pca_data <- pca_data[,apply(pca_data, 2, var, na.rm=TRUE) != 0]

# Compute the Principal Components
pca_data <- prcomp(pca_data, center = T, scale. = T)
summary(pca_data) # check variance

# Plotting PCA
library(ggbiplot)

#ggscreeplot(pca_data, type = "cev")

ggbiplot(pca_data,
         ellipse=F,
         choices=c(1,2),
         var.axes=F,
         labels=samples,
         groups=samples) +
  #scale_colour_manual(name="Origin", values= c("forest green", "red3", "dark blue"))+
  ggtitle("PCA of Mesoamerican Beans transcriptomics")+
  theme_minimal()+
  theme(legend.position = "bottom")


####################################
####################################
####################################

library(mixOmics)

neg_data_pls <- t(neg_data[,order(colnames(neg_data))])
#neg_data_pls <- neg_data_pls[,colnames(neg_data_pls) == "Module_13"]
tr_data_pls <- t(tr_data[,order(colnames(tr_data))])


liver.spls <- spls(tr_data_pls, neg_data_pls, ncomp = 2,
                   keepX = c(50), keepY= c(10), mode = "regression")

#liver.spls <- spls(tr_data_pls, neg_data_pls, ncomp = 2,
#                   keepX = c(50),  mode = "regression")

tune.spls <- perf(liver.spls, validation = "loo", progressBar = T)

tune.spls$measures$Q2.total$summary

plot(tune.spls$measures$Q2.total$summary$mean, ylim = c(-2,0.5),
     ylab = "Q2", xlab = "n of components")
abline(h = 0.0975)
abline(h = 0)


#list.keepX <- c(1:10, 15, 20, 30, 35 , 40 ,45, 50)
# tuning based on MAE
#tune.spls.MAE <- tune.spls(tr_data_pls, neg_data_pls, ncomp = 2,
#                           test.keepX = list.keepX,
#                           validation = "loo",
#                           progressBar = FALSE,
#                           measure = 'cor')
#plot(tune.spls.MAE, legend.position = 'topright')#

#tune.spls.MAE$choice.keepX


MySelectedVariables <- selectVar(liver.spls, comp = 1)
MySelectedVariables$Y$name # Selected genes on component 1
MySelectedVariables$X$name # Selected genes on component 1

MySelectedVariables$Y$name[order(MySelectedVariables$Y$name)]

sel_tr <- tr_data_raw[rownames(tr_data_raw) %in% MySelectedVariables$X$name,]
write.csv(sel_tr, "sel_tr.csv")

plotIndiv(liver.spls, comp = 1:2, rep.space= 'Y-variate', group = sub("_.*", "", rownames(neg_data_pls)),
          ind.names = sub("_.*", "", rownames(neg_data_pls)),
          legend = TRUE, title = 'Beans Wild vs Domesticate, sPLS comp 1 - 2, Y-space')

plotIndiv(liver.spls, comp = 1:2, rep.space= 'X-variate', group = sub("_.*", "", rownames(neg_data_pls)),
          ind.names = sub("_.*", "", rownames(neg_data_pls)),
          legend = TRUE, title = 'Beans Wild vs Domesticate, sPLS comp 1 - 2, X-space')

plotIndiv(liver.spls, comp = 1:2, rep.space= 'XY-variate', group = sub("_.*", "", rownames(neg_data_pls)),
          ind.names = sub("_.*", "", rownames(neg_data_pls)),
          legend = TRUE, title = 'Beans Wild vs Domesticate, sPLS comp 1 - 2, XY-space')

plotArrow(liver.spls, group=sub("_.*", "", rownames(neg_data_pls)), legend = TRUE,
          X.label = 'PLS comp 1', Y.label = 'PLS comp 2', legend.title = "Genepool")

plotVar(liver.spls, comp = 1:2)


# define red and green colors for the edges
color.edge <- color.GreenRed(50)  
# to save as a pdf
network(liver.spls, comp = 1:2, shape.node = c("rectangle", "rectangle"),
        color.node = c("white", "pink"), color.edge = color.edge, save = "pdf", name.save = "net")

cim(liver.spls, comp = 1, xlab = "Features", ylab = "Transcripts", save = "pdf", name.save = "cim",
    row.names = T, col.names = T)


