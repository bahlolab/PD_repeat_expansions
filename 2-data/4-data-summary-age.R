library(tidyr)
#load data

system("mkdir -p data_summary_counts/age/intermediate_files")
res = read.delim("../input_data/all_cohorts_FINAL_EH5_keep_loci_expanded_labelled_annotate_interruptions.txt")
cohort_size = read.delim("cohort_size/cohort_sizes_keep_only_by_age.txt")

cohort_size

#everything is already keep
#res = res[res$remove %in% "keep",] 


res= res[res$cohort %in% c("AMP-PD", "LBD"),]
locus_list = unique(res$locus)
res$age_group = ">=50"
res[res$age <50,]$age_group = "<50"


res_cohort = as.data.frame(table(res$locus, res$ancestry , res$cohort, res$cases, res$age_group, res$rep2_expanded))
colnames(res_cohort) = c("locus", "ancestry", "cohort", "cases", "age_group", "rep2_expanded", "count")

res_cohort = res_cohort[!res_cohort$rep2_expanded %in% "not_expanded",]
res_cohort_not = res_cohort
res_cohort_not$rep2_expanded = "not_expanded"

head(res_cohort_not)

res_cohort_not$uid = paste(res_cohort_not$ancestry, res_cohort_not$cohort, res_cohort_not$cases, res_cohort_not$age_group, sep= "_")
cohort_size$uid = paste(cohort_size$ancestry, cohort_size$cohort, cohort_size$cases, cohort_size$age_group, sep = "_")

res_cohort_not = merge(res_cohort_not, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort_not$count = res_cohort_not$count.y - res_cohort_not$count.x
res_cohort_not = res_cohort_not[,c(2:7, 10)]

res_cohort = rbind(res_cohort, res_cohort_not)
res_cohort_not = NULL

res_cohort$uid = paste(res_cohort$ancestry, res_cohort$cohort, res_cohort$cases, res_cohort$age_group, sep= "_")
res_cohort = merge(res_cohort, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort$percent = round( res_cohort$count.x/res_cohort$count.y*100, digits = 2)
res_cohort$uid = res_cohort$count.y = NULL
colnames(res_cohort)[7] = "count"
head(res_cohort) 




##############################

inter = res[res$has_intermediate %in% "yes",]
res_cohort_intermediate = as.data.frame(table(inter$locus, inter$ancestry , inter$cohort, inter$cases, inter$age_group, inter$rep2_expanded_intermediate))

colnames(res_cohort_intermediate) = c("locus", "ancestry", "cohort", "cases","age_group" ,"rep2_expanded_intermediate", "count")
res_cohort_intermediate = res_cohort_intermediate[!res_cohort_intermediate$rep2_expanded %in% "not_expanded",]
res_cohort_intermediate_not = res_cohort_intermediate

res_cohort_intermediate_not$rep2_expanded_intermediate = "not_expanded"

head(res_cohort_intermediate_not)

res_cohort_intermediate_not$uid = paste(res_cohort_intermediate_not$ancestry, res_cohort_intermediate_not$cohort, res_cohort_intermediate_not$cases, res_cohort_intermediate_not$age_group, sep= "_")
#cohort_size$uid = paste(cohort_size$ancestry, cohort_size$cohort, cohort_size$cases, sep = "_")

res_cohort_intermediate_not = merge(res_cohort_intermediate_not, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort_intermediate_not$count = res_cohort_intermediate_not$count.y - res_cohort_intermediate_not$count.x
res_cohort_intermediate_not = res_cohort_intermediate_not[,c(2:7, 10)]

res_cohort_intermediate = rbind(res_cohort_intermediate, res_cohort_intermediate_not)
res_cohort_intermediate_not = NULL

res_cohort_intermediate$uid = paste(res_cohort_intermediate$ancestry, res_cohort_intermediate$cohort, res_cohort_intermediate$cases, res_cohort_intermediate$age_group,  sep= "_")
res_cohort_intermediate = merge(res_cohort_intermediate, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort_intermediate$percent = round( res_cohort_intermediate$count.x/res_cohort_intermediate$count.y*100, digits = 2)
res_cohort_intermediate$uid = res_cohort_intermediate$count.y = NULL
colnames(res_cohort_intermediate)[7] = "count"
head(res_cohort_intermediate) 

##############################




write.table(res_cohort, file="data_summary_counts/age/intermediate_files/data_summary_all_cohorts_any_allele.txt", sep="\t", quote = F, row.names = F)
write.table(res_cohort_intermediate, file="data_summary_counts/age/intermediate_files/data_summary_all_cohorts_intermediate_any_allele.txt", sep="\t", quote = F, row.names = F)

############################################################
#recessive


res_cohort_recessive = as.data.frame(table(res$locus, res$ancestry , res$cohort, res$cases, res$age_group, res$both_expanded))
colnames(res_cohort_recessive) = c("locus", "ancestry", "cohort", "cases", "age_group", "both_expanded", "count")

res_cohort_recessive = res_cohort_recessive[!res_cohort_recessive$both_expanded %in% "no",]
res_cohort_recessive_not = res_cohort_recessive
res_cohort_recessive_not$both_expanded = "no"

head(res_cohort_recessive_not)

res_cohort_recessive_not$uid = paste(res_cohort_recessive_not$ancestry, res_cohort_recessive_not$cohort, res_cohort_recessive_not$cases, res_cohort_recessive_not$age_group, sep= "_")
#cohort_size$uid = paste(cohort_size$ancestry, cohort_size$cohort, cohort_size$cases, sep = "_")

res_cohort_recessive_not = merge(res_cohort_recessive_not, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort_recessive_not$count = res_cohort_recessive_not$count.y - res_cohort_recessive_not$count.x
res_cohort_recessive_not = res_cohort_recessive_not[,c(2:7, 10)]

res_cohort_recessive = rbind(res_cohort_recessive, res_cohort_recessive_not)
res_cohort_recessive_not = NULL


res_cohort_recessive$uid = paste(res_cohort_recessive$ancestry, res_cohort_recessive$cohort, res_cohort_recessive$cases, res_cohort_recessive$age_group, sep= "_")
res_cohort_recessive = merge(res_cohort_recessive, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort_recessive$percent = round( res_cohort_recessive$count.x/res_cohort_recessive$count.y*100, digits = 2)
res_cohort_recessive$uid = res_cohort_recessive$count.y = NULL
colnames(res_cohort_recessive)[7] = "count"
head(res_cohort_recessive) 

##############################





inter = res[res$has_intermediate %in% "yes",]
res_cohort_recessive_intermediate = as.data.frame(table(inter$locus, inter$ancestry , inter$cohort, inter$cases, inter$age_group, inter$both_expanded_intermediate))

colnames(res_cohort_recessive_intermediate) = c("locus", "ancestry", "cohort", "cases", "age_group", "both_expanded_intermediate", "count")
res_cohort_recessive_intermediate = res_cohort_recessive_intermediate[!res_cohort_recessive_intermediate$both_expanded %in% "no",]
res_cohort_recessive_intermediate_not = res_cohort_recessive_intermediate

res_cohort_recessive_intermediate_not$both_expanded_intermediate = "no"

head(res_cohort_recessive_intermediate_not)

res_cohort_recessive_intermediate_not$uid = paste(res_cohort_recessive_intermediate_not$ancestry, res_cohort_recessive_intermediate_not$cohort, res_cohort_recessive_intermediate_not$cases, res_cohort_recessive_intermediate_not$age_group, sep= "_")
#cohort_size$uid = paste(cohort_size$ancestry, cohort_size$cohort, cohort_size$cases, sep = "_")

res_cohort_recessive_intermediate_not = merge(res_cohort_recessive_intermediate_not, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort_recessive_intermediate_not$count = res_cohort_recessive_intermediate_not$count.y - res_cohort_recessive_intermediate_not$count.x
res_cohort_recessive_intermediate_not = res_cohort_recessive_intermediate_not[,c(2:7, 10)]

res_cohort_recessive_intermediate = rbind(res_cohort_recessive_intermediate, res_cohort_recessive_intermediate_not)
res_cohort_recessive_intermediate_not = NULL


res_cohort_recessive_intermediate$uid = paste(res_cohort_recessive_intermediate$ancestry, res_cohort_recessive_intermediate$cohort, res_cohort_recessive_intermediate$cases, res_cohort_recessive_intermediate$age_group, sep= "_")
res_cohort_recessive_intermediate = merge(res_cohort_recessive_intermediate, cohort_size[,c(5,6)], by = "uid", all.x = T)
res_cohort_recessive_intermediate$percent = round( res_cohort_recessive_intermediate$count.x/res_cohort_recessive_intermediate$count.y*100, digits = 2)
res_cohort_recessive_intermediate$uid = res_cohort_recessive_intermediate$count.y = NULL
colnames(res_cohort_recessive_intermediate)[7] = "count"
head(res_cohort_recessive_intermediate) 

##############################




write.table(res_cohort_recessive, file="data_summary_counts/age/intermediate_files/data_summary_all_cohorts_recessive.txt", sep="\t", quote = F, row.names = F)
write.table(res_cohort_recessive_intermediate, file="data_summary_counts/age/intermediate_files/data_summary_all_cohorts_intermediate_recessive.txt", sep="\t", quote = F, row.names = F)


##############################
#XLR

cohort_size = read.delim("cohort_size/cohort_sizes_keep_only_by_age_and_sex.txt")

res1 = res[res$inheritance %in% "XLR",]


res_cohort_xlr = as.data.frame(table(res1$locus, res1$ancestry, res1$cohort, res1$cases, res1$age_group, res1$sex, res1$both_expanded))
colnames(res_cohort_xlr) = c("locus", "ancestry", "cohort", "cases", "age_group", "sex", "both_expanded", "count")

res_cohort_xlr = res_cohort_xlr[!res_cohort_xlr$both_expanded %in% "no",]
res_cohort_xlr_not = res_cohort_xlr
res_cohort_xlr_not$both_expanded = "no"

head(res_cohort_xlr_not)

res_cohort_xlr_not$uid = paste(res_cohort_xlr_not$ancestry, res_cohort_xlr_not$cohort, res_cohort_xlr_not$sex, res_cohort_xlr_not$age_group, res_cohort_xlr_not$cases, sep= "_")
cohort_size$uid = paste(cohort_size$ancestry, cohort_size$cohort, cohort_size$sex,  cohort_size$age_group, cohort_size$cases, sep = "_")

res_cohort_xlr_not = merge(res_cohort_xlr_not, cohort_size[,c(6,7)], by = "uid", all.x = T)
res_cohort_xlr_not$count = res_cohort_xlr_not$count.y - res_cohort_xlr_not$count.x
res_cohort_xlr_not = res_cohort_xlr_not[,c(2:8, 11)]

res_cohort_xlr = rbind(res_cohort_xlr, res_cohort_xlr_not)
res_cohort_xlr_not = NULL


res_cohort_xlr$uid = paste(res_cohort_xlr$ancestry, res_cohort_xlr$cohort, res_cohort_xlr$sex,  res_cohort_xlr$age_group,  res_cohort_xlr$cases, sep= "_")
res_cohort_xlr = merge(res_cohort_xlr, cohort_size[,c(6,7)], by = "uid", all.x = T)
res_cohort_xlr$percent = round( res_cohort_xlr$count.x/res_cohort_xlr$count.y*100, digits = 2)
res_cohort_xlr$uid = res_cohort_xlr$count.y = NULL
colnames(res_cohort_xlr)[7] = "count"
head(res_cohort_xlr) 




write.table(res_cohort_xlr, file="data_summary_counts/age/data_summary_all_cohorts_XLR_by_sex.txt", sep="\t", quote = F, row.names = F)

#note, XLR does not have intermediate range

############ merge data




res1 = read.delim("data_summary_counts/age/intermediate_files/data_summary_all_cohorts_any_allele.txt")
res2 = read.delim("data_summary_counts/age/intermediate_files/data_summary_all_cohorts_intermediate_any_allele.txt")
res3 = read.delim("data_summary_counts/age/intermediate_files/data_summary_all_cohorts_recessive.txt")
res4 = read.delim("data_summary_counts/age/intermediate_files/data_summary_all_cohorts_intermediate_recessive.txt")

res1$rep2_expanded = gsub("not_expanded", "no", res1$rep2_expanded)
res1$rep2_expanded = gsub("expanded", "yes", res1$rep2_expanded)
res2$rep2_expanded_intermediate = gsub("not_expanded", "no", res2$rep2_expanded_intermediate)
res2$rep2_expanded_intermediate = gsub("expanded", "yes", res2$rep2_expanded_intermediate)

res1$uid = paste(res1$locus, res1$ancestry, res1$cohort, res1$cases, res1$age_group, res1$rep2_expanded, sep = "_")
res2$uid = paste(res2$locus, res2$ancestry, res2$cohort, res2$cases, res2$age_group, res2$rep2_expanded_intermediate, sep = "_")
res3$uid = paste(res3$locus, res3$ancestry, res3$cohort, res3$cases, res3$age_group, res3$both_expanded, sep = "_")
res4$uid = paste(res4$locus, res4$ancestry, res4$cohort, res4$cases, res4$age_group, res4$both_expanded_intermediate, sep = "_")

colnames(res1)[c(7,8)] = paste0(colnames(res1)[c(7,8)], "_DOM" )
colnames(res2)[c(7,8)] = paste0(colnames(res2)[c(7,8)], "_DOMi" )
colnames(res3)[c(7,8)] = paste0(colnames(res3)[c(7,8)], "_REC" )
colnames(res4)[c(7,8)] = paste0(colnames(res4)[c(7,8)], "_RECi" )

res = merge(res1, res2[,c(7,8,9)], by="uid", all = T)
res = merge(res, res3[,c(7,8,9)], by="uid", all = T)
res = merge(res, res4[,c(7,8,9)], by="uid", all = T)

res$uid = NULL
colnames(res)[6] = "expanded"
write.table(res, file="data_summary_counts/age/data_summary_all_cohort.txt", sep="\t", quote = F, row.names = F)











