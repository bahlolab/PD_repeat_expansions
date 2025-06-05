system("mkdir cohort_size")

res = read.delim("../input_data/all_cohorts_FINAL_EH5_keep_loci_expanded_labelled_annotate_interruptions.txt")

uniq = res[!duplicated(res$sample),]
uniq = uniq[!(uniq$cohort %in% "gnomad" & uniq$gene != "TCF4"),]
cohort_all_ethnicities    = as.data.frame(table( uniq$cohort, uniq$ancestry, uniq$cases))
subcohort_all_ethnicities = as.data.frame(table( uniq$subcohort, uniq$ancestry, uniq$cases))

colnames(cohort_all_ethnicities) = c("cohort", "ancestry", "cases", "count")
colnames(subcohort_all_ethnicities) = c("cohort", "ancestry", "cases", "count")

write.table(cohort_all_ethnicities, file="cohort_size/cohort_sizes_keep_only.txt", sep="\t", row.names = F, quote = F)
write.table(subcohort_all_ethnicities, file="cohort_size/subcohort_sizes_keep_only.txt", sep="\t", row.names = F, quote = F)



#split by sex


cohort_all_ethnicities    = as.data.frame(table( uniq$cohort, uniq$ancestry, uniq$sex, uniq$cases))
subcohort_all_ethnicities = as.data.frame(table( uniq$subcohort, uniq$ancestry, uniq$sex, uniq$cases))

colnames(cohort_all_ethnicities) = c("cohort", "ancestry", "sex", "cases", "count")
colnames(subcohort_all_ethnicities) = c("cohort", "ancestry", "sex","cases", "count")

write.table(cohort_all_ethnicities, file="cohort_size/cohort_sizes_keep_only_by_sex.txt", sep="\t", row.names = F, quote = F)
write.table(subcohort_all_ethnicities, file="cohort_size/subcohort_sizes_keep_only_by_sex.txt", sep="\t", row.names = F, quote = F)


#split by age

head(uniq)
uniq$age_group = ">=50"
uniq = uniq[uniq$cohort %in% c("AMP-PD", "LBD"),]
uniq[uniq$age < 50,]$age_group = "<50"

cohort_all_ethnicities    = as.data.frame(table( uniq$cohort, uniq$ancestry, uniq$age_group, uniq$cases))
subcohort_all_ethnicities = as.data.frame(table( uniq$subcohort, uniq$ancestry, uniq$age_group, uniq$cases))

colnames(cohort_all_ethnicities) = c("cohort", "ancestry", "age_group", "cases", "count")
colnames(subcohort_all_ethnicities) = c("cohort", "ancestry", "age_group","cases", "count")

write.table(cohort_all_ethnicities, file="cohort_size/cohort_sizes_keep_only_by_age.txt", sep="\t", row.names = F, quote = F)
write.table(subcohort_all_ethnicities, file="cohort_size/subcohort_sizes_keep_only_by_age.txt", sep="\t", row.names = F, quote = F)

#split by age and sex



cohort_all_ethnicities    = as.data.frame(table( uniq$cohort, uniq$ancestry, uniq$age_group, uniq$sex,  uniq$cases))
subcohort_all_ethnicities = as.data.frame(table( uniq$subcohort, uniq$ancestry, uniq$age_group, uniq$sex,  uniq$cases))

colnames(cohort_all_ethnicities) = c("cohort", "ancestry", "age_group", "sex", "cases", "count")
colnames(subcohort_all_ethnicities) = c("cohort", "ancestry", "age_group", "sex", "cases", "count")

write.table(cohort_all_ethnicities, file="cohort_size/cohort_sizes_keep_only_by_age_and_sex.txt", sep="\t", row.names = F, quote = F)
write.table(subcohort_all_ethnicities, file="cohort_size/subcohort_sizes_keep_only_by_age_and_sex.txt", sep="\t", row.names = F, quote = F)

