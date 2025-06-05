meta = read.delim("../input_data/all_metadata_expanded_v2_clinical.txt")


res0 = read.delim("../input_data/all_cohorts_FINAL_EH5_keep_loci_expanded_labelled_annotate_interruptions.txt")
res = res0[res0$gene %in% "ATXN2",]
res = res[res$cohort %in% c("AMP-PD", "LBD"),]
res = res[res$ancestry %in% "EUR",]
res = res[res$rep2 > 32,]
res = res[,c(1,6,7,8,9,10,13, 14,16, 25, 26,29,30)]

meta = meta[meta$sample %in% res$sample,]

meta = meta[,c(1,5, 12:14,16:19, 29:33)]

res = merge(res, meta, by="sample")

res = res[ order(-res$rep2),]
res$ID = paste0("Patient ", 1:nrow(res))

res = as.data.frame( t(res) )

res$measurments= row.names(res)

res = res[,c(7, 1:6)]
row.names(res) = 1:nrow(res)

write.table(res, file="Table_1_ATXN2.txt", sep="\t", quote = F, row.names = F, col.names = F)

res = read.delim("Table_1_ATXN2.txt")



