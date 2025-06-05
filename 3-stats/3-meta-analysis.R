library(poolr)

res = read.delim("fishers_test_EUR.txt")


fres = NULL
mres = NULL
for(inheritance in unique(res$inheritance)){
for(locus in unique(res$locus)){
  test = res[res$locus %in% locus & res$inheritance %in% inheritance & res$n2 != 0 ,]
  
  x = fisher(test$p)
  fres = rbind(fres,  c(locus, inheritance, round(x$p, digits = 3), x$k))
  
  mx = matrix(colSums(test[,9:12]), nrow = 2)
  mxf = fisher.test(mx)
  mres = rbind(mres,  c(locus, inheritance, mxf$estimate, round(mxf$p.value, digits = 3), 
                        mxf$conf.int[1], mxf$conf.int[2], mx[1], mx[2], mx[3], mx[4]))
  
  }}

fres = as.data.frame(fres)
colnames(fres) = c("locus", "inheritance", "p", "number")
fres = fres[fres$number > 0,]

mres = as.data.frame(mres)
colnames(mres) = colnames(res)[c(2,4:12)]
mres = mres[mres$OR > 0,]

