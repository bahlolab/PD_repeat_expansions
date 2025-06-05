library(tidyr)
library(dplyr)
#load data

system("mkdir -p data_summary_counts/intermediate_files")
res = read.delim("../input_data/all_cohorts_FINAL_EH5_keep_loci_expanded_labelled_annotate_interruptions.txt")
cohort_size = read.delim("cohort_size/cohort_sizes_keep_only.txt")
subcohort_size = read.delim("cohort_size/subcohort_sizes_keep_only.txt")

cohort_size = cohort_size[cohort_size$ancestry %in% "EUR",]
table(res[res$cohort %in% "AMP-PD" & !duplicated(res$sample),]$ancestry)
#everything is already keep
#res = res[res$remove %in% "keep",] 

####EUR only
####
res = res[res$ancestry %in% "EUR",]

#exclude benign or not relevant loci
exclude = c("SCA31_BEAN1_AAAAT",  "FECD3_TCF4")
res = res[!res$locus %in% exclude, ]

locus_list = unique(res$locus)

#pathogenic carrier, dominant or 2 copies recessive

res$pathogenic = "no"
unique(res$inheritance)
res[res$inheritance %in% c("AD", "XLD") & res$rep2_expanded %in% "expanded",]$pathogenic = "yes"
res[res$inheritance %in% c("AR", "XLR") & res$both_expanded %in% "yes",]$pathogenic = "yes"

ataxia = res[res$locus %in% 
               c("SCA1_ATXN1", "SCA17_TBP", "SCA2_ATXN2", "SCA4_ZFHX3",
                 "SCA8_ATXN8OS", "SCA27B_FGF14",
                 "SCA10_ATXN10", "SCA12_PPP2R2B", "SCA3_ATXN3", "FRDA_FXN",
                 "SCA36_NOP56", "SCA6_CACNA1A", "SCA7_ATXN7", "CANVAS_RFC1_AAGGG", "CANVAS_RFC1_ACAGG"),]

res_cohort = as.data.frame(table(ataxia$sample , ataxia$cohort, ataxia$cases, ataxia$pathogenic))
colnames(res_cohort) = c("sample",  "cohort", "cases", "pathogenic", "count")

res_cohort = res_cohort[!res_cohort$pathogenic %in% "no",]
res_cohort = res_cohort[!res_cohort$sample %in% "",]

res_cohort = res_cohort[res_cohort$count > 0 ,]

#this is the number of people with at least one expansion
res_cohort = res_cohort %>% group_by(cohort, cases, pathogenic) %>%
  summarise(count = n()) 

#number of people with no expansions
res_cohort_not = res_cohort
res_cohort_not$pathogenic = "no"

head(res_cohort_not)

res_cohort_not$uid = paste(res_cohort_not$cohort, res_cohort_not$cases, sep= "_")
cohort_size$uid = paste(cohort_size$cohort, cohort_size$cases, sep = "_")

res_cohort_not = merge(res_cohort_not, cohort_size[,c(4,5)], by = "uid", all.x = T)
res_cohort_not$count = res_cohort_not$count.y - res_cohort_not$count.x
res_cohort_not = res_cohort_not[,c(2:4, 7)]

res_cohort = rbind(res_cohort, res_cohort_not)
res_cohort_not = NULL

res_cohort$uid = paste(res_cohort$cohort, res_cohort$cases, sep= "_")
res_cohort = merge(res_cohort, cohort_size[,c(4,5)], by = "uid", all.x = T)
res_cohort$percent = round( res_cohort$count.x/res_cohort$count.y*100, digits = 2)
res_cohort$uid = res_cohort$count.y = NULL
colnames(res_cohort)[4] = "count"
head(res_cohort) 







not_ataxia = res[!res$locus %in% 
                   c("SCA1_ATXN1", "SCA17_TBP", "SCA2_ATXN2", "SCA4_ZFHX3",
                     "SCA8_ATXN8OS", "SCA27B_FGF14",
                     "SCA10_ATXN10", "SCA12_PPP2R2B", "SCA3_ATXN3", 
                     "SCA36_NOP56", "SCA6_CACNA1A", "SCA7_ATXN7", "CANVAS_RFC1_AAGGG", "CANVAS_RFC1_ACAGG"),]

res_cohort_nonataxia = as.data.frame(table(not_ataxia$sample , not_ataxia$cohort, not_ataxia$cases, not_ataxia$pathogenic))
colnames(res_cohort_nonataxia) = c("sample",  "cohort", "cases", "pathogenic", "count")

res_cohort_nonataxia = res_cohort_nonataxia[!res_cohort_nonataxia$pathogenic %in% "no",]
res_cohort_nonataxia = res_cohort_nonataxia[!res_cohort_nonataxia$sample %in% "",]

res_cohort_nonataxia = res_cohort_nonataxia[res_cohort_nonataxia$count > 0 ,]

#this is the number of people with at least one expansion
res_cohort_nonataxia = res_cohort_nonataxia %>% group_by(cohort, cases, pathogenic) %>%
  summarise(count = n()) 

#number of people with no expansions
res_cohort_nonataxia_not = res_cohort_nonataxia
res_cohort_nonataxia_not$pathogenic = "no"

head(res_cohort_nonataxia_not)

res_cohort_nonataxia_not$uid = paste(res_cohort_nonataxia_not$cohort, res_cohort_nonataxia_not$cases, sep= "_")
cohort_size$uid = paste(cohort_size$cohort, cohort_size$cases, sep = "_")

res_cohort_nonataxia_not = merge(res_cohort_nonataxia_not, cohort_size[,c(4,5)], by = "uid", all.x = T)
res_cohort_nonataxia_not$count = res_cohort_nonataxia_not$count.y - res_cohort_nonataxia_not$count.x
res_cohort_nonataxia_not = res_cohort_nonataxia_not[,c(2:4, 7)]

res_cohort_nonataxia = rbind(res_cohort_nonataxia, res_cohort_nonataxia_not)
res_cohort_nonataxia_not = NULL

res_cohort_nonataxia$uid = paste(res_cohort_nonataxia$cohort, res_cohort_nonataxia$cases, sep= "_")
res_cohort_nonataxia = merge(res_cohort_nonataxia, cohort_size[,c(4,5)], by = "uid", all.x = T)
res_cohort_nonataxia$percent = round( res_cohort_nonataxia$count.x/res_cohort_nonataxia$count.y*100, digits = 2)
res_cohort_nonataxia$uid = res_cohort_nonataxia$count.y = NULL
colnames(res_cohort_nonataxia)[4] = "count"
head(res_cohort_nonataxia) 


res_cohort_all_loci = as.data.frame(table(res$sample , res$cohort, res$cases, res$pathogenic))
colnames(res_cohort_all_loci) = c("sample",  "cohort", "cases", "pathogenic", "count")

res_cohort_all_loci = res_cohort_all_loci[!res_cohort_all_loci$pathogenic %in% "no",]
res_cohort_all_loci = res_cohort_all_loci[!res_cohort_all_loci$sample %in% "",]

res_cohort_all_loci = res_cohort_all_loci[res_cohort_all_loci$count > 0 ,]

#this is the number of people with at least one expansion
res_cohort_all_loci = res_cohort_all_loci %>% group_by(cohort, cases, pathogenic) %>%
  summarise(count = n()) 

#number of people with no expansions
res_cohort_all_loci_not = res_cohort_all_loci
res_cohort_all_loci_not$pathogenic = "no"

head(res_cohort_all_loci_not)

res_cohort_all_loci_not$uid = paste(res_cohort_all_loci_not$cohort, res_cohort_all_loci_not$cases, sep= "_")
cohort_size$uid = paste(cohort_size$cohort, cohort_size$cases, sep = "_")

res_cohort_all_loci_not = merge(res_cohort_all_loci_not, cohort_size[,c(4,5)], by = "uid", all.x = T)
res_cohort_all_loci_not$count = res_cohort_all_loci_not$count.y - res_cohort_all_loci_not$count.x
res_cohort_all_loci_not = res_cohort_all_loci_not[,c(2:4, 7)]

res_cohort_all_loci = rbind(res_cohort_all_loci, res_cohort_all_loci_not)
res_cohort_all_loci_not = NULL

res_cohort_all_loci$uid = paste(res_cohort_all_loci$cohort, res_cohort_all_loci$cases, sep= "_")
res_cohort_all_loci = merge(res_cohort_all_loci, cohort_size[,c(4,5)], by = "uid", all.x = T)
res_cohort_all_loci$percent = round( res_cohort_all_loci$count.x/res_cohort_all_loci$count.y*100, digits = 2)
res_cohort_all_loci$uid = res_cohort_all_loci$count.y = NULL
colnames(res_cohort_all_loci)[4] = "count"
head(res_cohort_all_loci) 


res_cohort_all_loci$version = "all_loci"
res_cohort$version = "ataxia"
res_cohort_nonataxia$version = "non-ataxia"

results = rbind(res_cohort, res_cohort_all_loci, res_cohort_nonataxia)

write.table(results, file="data_summary_counts/data_summary_all_merged_EUR.txt", sep= "\t", quote = F, row.names = F)




