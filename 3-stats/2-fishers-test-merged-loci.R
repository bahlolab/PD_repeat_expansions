library(ggplot2)
#system("mkdir plots")

x = read.delim("../data/data_summary_counts/data_summary_all_merged_EUR.txt")

head(x)


cohorts = unique(x$cohort)
cohorts = cohorts[cohorts != "gnomad"]
loci = unique(x$version)


res = NULL

# Cases
for(case in c("Case", "Other"))
    for(cohort in cohorts) {
      for(locus in loci) {
        
        x1 = x[x$cohort %in% cohort & x$version %in% locus,]
        if(nrow(x1[x1$cases %in% case,]) == 2){
          n1 = x1[x1$cases %in% case & x1$pathogenic %in% "yes",]$count
          n2 = x1[x1$cases %in% case & x1$pathogenic %in% "no",]$count
          n3 = x1[x1$cases %in% "Control" & x1$pathogenic %in% "yes",]$count
          n4 = x1[x1$cases %in% "Control" & x1$pathogenic %in% "no",]$count
          mx = c(n1, n2, n3, n4)
          mx = matrix(mx, ncol = 2)
          
          mx_f = fisher.test(mx)
          res = rbind(res, c(cohort, locus, paste0(case, "-control") , c(mx_f$estimate, mx_f$p.value, mx_f$conf.int[1], mx_f$conf.int[2] ), n1, n2, n3, n4 ))
        }else{}
        }}




res = as.data.frame(res)
colnames(res) = c("cohort", "locus", "comparison", "OR", "p", "CIL", "CIU", "n1", "n2", "n3", "n4")


write.table(res, file="fishers_test_merged_loci_EUR.txt", sep = "\t", quote = F, row.names = F)
