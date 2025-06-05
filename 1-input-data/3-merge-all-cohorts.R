#Merge genotypes from AMP-PD and gnomAD
amp = read.delim("all_results_AMPPD_EH5_keep_loci_expanded_labelled.txt")
gnomad = read.delim("all_results_gnomad_EH5_keep_loci_expanded_labelled.txt")

amp = amp[,-c(5:14)]
amp$cohort2 = "AMP-PD"

amp[amp$cohort %in% "LBD",]$cohort2 = "LBD"

colnames(amp)[11] = "age"



head(gnomad)

setdiff(colnames(amp), colnames(gnomad))

gnomad$cases = "Control"
gnomad$history = "unknown"
gnomad$group = gnomad$cases
gnomad$remove = "keep"
gnomad$race = NA
gnomad$cohort2 = gnomad$cohort


gnomad = gnomad[,c("sample",  "locus", "rep1", "rep2", "cohort", "cases", "history", "group", 
                   "remove", "sex", "age" ,  "race" ,   "population",
                   "rep1_expanded", "rep2_expanded", "has_intermediate", "rep1_expanded_intermediate", "rep2_expanded_intermediate", "cohort2" )]



colnames(gnomad)[13] = "ancestry"
gnomad$age = NA
gnomad$ancestry = toupper(gnomad$ancestry)
gnomad$ancestry = gsub("NFE", "EUR", gnomad$ancestry)

res = rbind(amp, gnomad)
  

#add both_expanded for recessive hits

res$both_expanded = "no"
res$both_expanded_intermediate = "no"
res[res$rep1_expanded %in% "expanded" & res$rep2_expanded %in% "expanded", ]$both_expanded = "yes"
res[res$rep1_expanded_intermediate %in% "expanded" & res$rep2_expanded_intermediate %in% "expanded", ]$both_expanded_intermediate = "yes"


res = res[,c(1:5, 19, 6:18, 20:21)]
colnames(res)[5] = "subcohort"
colnames(res)[6] = "cohort"

res$locus = gsub("FGF14_AAG", "SCA27B_FGF14", res$locus)
res$gene = gsub(".*_", "", res$locus)
res[res$locus %in% "THAP11_AGC",]$gene = "THAP11"
res[res$locus %in% "CANVAS_RFC1_AAGGG",]$gene = "RFC1_AAGGG"
res[res$locus %in% "CANVAS_RFC1_ACAGG",]$gene = "RFC1_ACAGG"
res$disorder = gsub("_.*", "", res$locus)
res = res[,c(1,2,22,23,3:21)]

db = read.delim("exSTRa_hg38.txt", skip=1)
db$locus = gsub("FGF14_AAG", "SCA27B_FGF14", db$locus)

res = merge(res, db[,c(1,4)], by= "locus", all.x=T)

res = res[,c(2,1,3:24)]

#For males, change both expanded to YES if they are expanded on one Y allele - for XL only.

head(res)

res[res$sex %in% "male" & res$inheritance %in% c("XLD", "XLR") & res$rep2_expanded %in% "expanded",]$both_expanded = "yes"
res[res$sex %in% "male" & res$inheritance %in% c("XLD", "XLR") & res$rep2_expanded_intermediate %in% "expanded",]$both_expanded_intermediate = "yes"

#Exclude loci not relevant to the study
removed = res[!(res$remove %in% "keep"),]
res = res[res$remove %in% "keep",]

write.table(res, file="all_cohorts_FINAL_EH5_keep_loci_expanded_labelled.txt", sep="\t", row.names = F, quote = F)
write.table(removed, file="all_results_AMPPD_EH5_keep_loci_expanded_labelled_REMOVED_ONLY.txt", sep="\t", row.names = F, quote = F)


