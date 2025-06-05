library(ggplot2)
#system("mkdir plots")

x = read.delim("../data/data_summary_counts/data_summary_all_cohort.txt")

head(x)

x = x[x$ancestry %in% "EUR" ,]
x = x[,c(1,3:6,8,10,12)] 
head(x)

cohorts = unique(x$cohort)
cohorts = cohorts[cohorts != "gnomad"]
loci = unique(x$locus)
loci_int = unique(x[ !is.na(x$count_DOMi), ]$locus)

res = NULL

# Cases
for(case in c("Case", "Other"))
  for(model in c(5:8)){
    for(cohort in cohorts) {
      for(locus in loci) {
        
        x1 = x[x$cohort %in% cohort & x$locus %in% locus,]
        x1 = x1[ !is.na(x1[model]), ]
        
        if(nrow(x1) == 0){} else{
          
          n1 = unlist(unname(x1[x1$cases %in% case & x1$expanded %in% "yes",][model]))
          n2 = unlist(unname(x1[x1$cases %in% case & x1$expanded %in% "no",][model]))
          n3 = unlist(unname(x1[x1$cases %in% "Control" & x1$expanded %in% "yes",][model]))
          n4 = unlist(unname(x1[x1$cases %in% "Control" & x1$expanded %in% "no",][model]))
          mx = c(n1, n2, n3, n4)
          mx = matrix(mx, ncol = 2)
          
          mx_f = fisher.test(mx)
          res = rbind(res, c(cohort, locus, paste0(case, "-control") ,colnames(x1[model]), c(mx_f$estimate, mx_f$p.value, mx_f$conf.int[1], mx_f$conf.int[2] ), n1, n2, n3, n4 ))
        }}}}




res = as.data.frame(res)
colnames(res) = c("cohort", "locus", "comparison", "inheritance", "OR", "p", "CIL", "CIU", "n1", "n2", "n3", "n4")

res$inheritance = gsub("count_", "", res$inheritance)
res$inheritance = gsub("DOMi", "dominant (intermediate)", res$inheritance)
res$inheritance = gsub("DOM", "dominant", res$inheritance)
res$inheritance = gsub("RECi", "recessive (intermediate)", res$inheritance)
res$inheritance = gsub("REC", "recessive", res$inheritance)

write.table(res, file="fishers_test_EUR.txt", sep = "\t", quote = F, row.names = F)


