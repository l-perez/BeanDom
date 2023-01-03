

data <- as.matrix(read.csv("00_DATA_normalized_PQN.csv", row.names = 1))
cor_data <- cor(data)

temp <- unlist(strsplit(rownames(cor_data), "mz."))
mz_data <- temp[seq(1, length(temp), 2)]
rt_data <- temp[which(!temp %in% mz_data)]
rm(temp)

mz_data <- as.numeric(sub(".*?\\.", "", mz_data))
rt_data <- as.numeric(sub("min", "", rt_data))

diff_rt <-abs(outer(rt_data, rt_data, FUN = "-"))
rownames(diff_rt) <- rownames(cor_data)
colnames(diff_rt) <- colnames(cor_data)

group_table <- data.frame(feature = rownames(cor_data), group = NA)

for(n in group_table$feature[1]){
  rt_select <- diff_rt[which(rownames(diff_rt) %in% n),]
  rt_select <- rt_select[which(rt_select < 0.06)]
  corr_select <- cor_data[which(rownames(cor_data) %in% n),]
  corr_select <- corr_select[which(corr_select > 0.8)]
  group_select <- names(rt_select)[which(names(rt_select) %in% names(corr_select))]
  group_table$group[which(group_table$feature %in% group_select)] <- 1
}

for(n in group_table$feature){
  if(is.na(group_table$group[which(group_table$feature %in% n)])){
    rt_select <- diff_rt[which(rownames(diff_rt) %in% n),]
    rt_select <- rt_select[which(rt_select < 0.06)]
    corr_select <- cor_data[which(rownames(cor_data) %in% n),]
    corr_select <- corr_select[which(corr_select > 0.8)]
    group_select <- names(rt_select)[which(names(rt_select) %in% names(corr_select))]
    if(all(is.na(group_table$group[which(group_table$feature %in% group_select)]))){
      group_table$group[which(group_table$feature %in% group_select)] <- max(group_table$group, na.rm = T)+1
    } else {group_table$group[which(group_table$feature %in% group_select)] <- max(group_table$group[which(group_table$feature %in% group_select)], na.rm = T) }
    
  }
}

write.csv(group_table, "feature_group.csv", row.names = F)

group_table$median <- apply(data, 2, median)

library(dplyr)
selected_vars <- group_table %>% group_by(group) %>% top_n(1, median)
write.csv(selected_vars, "sel_median.csv")

selected_data <- data[, which(colnames(data) %in% selected_vars$feature)]
write.csv(selected_data, "00_DATA_GroupedData.csv")




