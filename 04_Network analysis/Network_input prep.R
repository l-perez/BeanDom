
# Load session from 03_Data analysis
load("210715_nets.Rdata")


# fix names cor
cornet_neg$Feature_A <- simnet_neg$Feature_A
cornet_neg$Feature_B <- simnet_neg$Feature_B

cornet_pos$Feature_A <- simnet_pos$Feature_A
cornet_pos$Feature_B <- simnet_pos$Feature_B

# filter by correlation and similarity
cornet_neg <- cornet_neg[cornet_neg$spearman_cor > 0.65,]
cornet_pos <- cornet_pos[cornet_pos$spearman_cor > 0.65,]

simnet_neg <- simnet_neg[simnet_neg$Similarity_score > 0.65,]
simnet_pos <- simnet_pos[simnet_pos$Similarity_score > 0.65,]


#write to .csv
write.csv(bionet_neg, "net_bio_neg.csv", row.names = F)
write.csv(bionet_pos, "net_bio_pos.csv", row.names = F)

write.csv(cornet_neg, "net_cor_neg.csv", row.names = F)
write.csv(cornet_pos, "net_cor_pos.csv", row.names = F)

write.csv(simnet_neg, "net_sim_neg.csv", row.names = F)
write.csv(simnet_pos, "net_sim_pos.csv", row.names = F)

# select and mark features previously identified in other work
temp <- neg_feat[,c(1, 27)]
temp_2 <- neg_feat$Class_tpj
temp_2[is.na(temp_2)] <- "Unknown"
temp$class <- temp_2

write.csv(temp, "features_neg.csv", row.names = F)

temp <- pos_feat[,c(1, 27)]
temp_2 <- pos_feat$Class_tpj
temp_2[is.na(temp_2)] <- "Unknown"
temp$class <- temp_2

write.csv(temp, "features_pos.csv", row.names = F)






