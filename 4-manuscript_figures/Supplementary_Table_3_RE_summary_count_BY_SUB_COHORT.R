library(ggplot2)
#system("mkdir plots")

x = read.delim("../data/data_summary_counts/data_summary_all_subcohort.txt")
res = read.delim("../input_data/all_cohorts_FINAL_EH5_keep_loci_expanded_labelled.txt")

map = res[,c("locus", "gene", "disorder")]
map = map[ !duplicated(map$locus),]

head(x)

x = x[x$ancestry %in% "EUR" ,]

cohort_size = read.delim("../data/cohort_size/subcohort_sizes_keep_only.txt")
cohort_size = cohort_size[cohort_size$ancestry %in% "EUR",]

cohort_size$uid = paste0(cohort_size$cohort, "-", cohort_size$cases)
x$uid = paste0(x$subcohort, "-", x$cases)

x = merge(x, cohort_size[,c(4:5)], by = "uid")
x = x[x$expanded %in% "yes",]

x$ancestry = x$expanded = x$uid = NULL


x$dom = paste0(x$count_DOM, " (", x$percent_DOM, "%)")
x$domi = paste0(x$count_DOMi, " (", x$percent_DOMi, "%)")
x$rec = paste0(x$count_REC, " (", x$percent_REC, "%)")
x$reci = paste0(x$count_RECi, " (", x$percent_RECi, "%)")


x[grepl("NA.*NA", x$reci),]$reci = "-"
x[grepl("NA.*NA", x$domi),]$domi = "-"



x = x[,c(1:3, 12:16)]


x = x[!x$subcohort %in% "UKBB-PD",]
x = x[!(x$subcohort %in% "gnomad" & x$cases %in% "Case"),] ; head(x)
x = x[!(x$subcohort %in% "gnomad" & x$cases %in% "Other"),] ; head(x)
x = x[!(x$subcohort %in% "LBD" & x$cases %in% "Case"),] ; head(x)

x$subcohort = gsub("gnomad", "gnomAD", x$subcohort)

x = merge(map, x, by= "locus")

x = x[order(x$locus, x$subcohort, x$cases),]
x$locus = NULL
colnames(x) = c("gene", "disorder", "cohort", "cases", "count", 
                "RE on longer allele", "RE on longer allele (inter)", 
                "RE on both alleles", "RE on both alleles (inter)" )

write.table(x, file="Supplementary_Table_3_RE_summary_count_SUB_COHORT.txt", sep="\t", quote = F, row.names = F, col.names = T)


