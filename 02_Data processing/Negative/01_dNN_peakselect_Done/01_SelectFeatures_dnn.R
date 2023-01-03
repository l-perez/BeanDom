

# input data ----------------------------------------------------------------------
all_content = readLines("00_DATA_bean_metaboanalyst_genpool_MZmineNofilter.csv")
data = read.csv(textConnection(all_content[-2]), row.names = 1)

select_feature <- read.csv("00_Subset_dNN_WindowMaxIntensityTable(20200414h1020)_subset(20200529h1714).csv")[1:3]

temp <- unlist(strsplit(rownames(data), "mz."))
mz_data <- temp[seq(1, length(temp), 2)]
rt_data <- temp[which(!temp %in% mz_data)]
rm(temp)

mz_data <- as.numeric(sub(".*/", "", mz_data))
rt_data <- as.numeric(sub("min", "", rt_data))

diff_mz <-as.data.frame(abs(outer(mz_data, select_feature$windMzMid, FUN = "-")))
diff_rt <-as.data.frame(abs(outer(rt_data, select_feature$windRtMid, FUN = "-")))

diff_mz <- as.data.frame(diff_mz < 0.5)
diff_rt <- as.data.frame(diff_rt < 0.05)

feature_match <- diff_mz*diff_rt

clean_data <- data[as.logical(rowSums(feature_match != 0)), ]

colnames(clean_data) <- sub("X", "", colnames(clean_data))
write.csv(clean_data, "00_DATA_bean_metaboanalyst_filterdNN.csv")
