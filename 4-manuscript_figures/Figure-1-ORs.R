library(ggplot2)
library(dplyr)
library(tidyr)

df = read.delim(file="../stats/fishers_test_EUR.txt")

df = df[!df$locus %in% "FECD3_TCF4",]


df = df[!(df$cohort %in% "LBD" & df$comparison %in% "Case-control"),]
df = df[!(df$cohort %in% "UKBB-PD" ),]

df = df[!(df$n1 == 0 & df$n3 == 0) , ]

df$fdr = p.adjust(df$p, method = "fdr")

gsub("", "", df$comparison)

df$logOR = log2(df$OR)
df$logCIU = log2(df$CIU)
df$logCIL = log2(df$CIL)
df = df[order(df$locus),]
df$locus = factor(df$locus, levels = rev(unique(df$locus)))
df$fdr_p = "> 0.05"
#df[df$fdr < 0.05,]$fdr_p = "< 0.05" 

df$threshold = "full RE"
df[grepl("inter", df$inheritance),]$threshold = "intermediate RE"
df$threshold

df$inheritance =gsub(" [(]intermediate[)]", "", df$inheritance)

df = df %>% separate_wider_delim(locus, "_", names = c("disorder", "gene", "motif"), too_few = "align_start" )

df$locus = df$gene

df[ !is.na(df$motif) ,]$locus = paste0(df[ !is.na(df$motif) ,]$gene, " (", df[ !is.na(df$motif) ,]$motif, ")")

df = df[order(df$locus),]
df$locus = factor(df$locus, levels = (unique(df$locus)))

db = read.delim("../input_data/exSTRa_hg38.txt", skip = 1)
db = db[db$inheritance %in% c("AR",  "XLR"),]
db = db[,c(1,4)]

db = db %>% separate_wider_delim(locus, "_", names = c("disorder", "gene", "motif"), too_few = "align_start" )

db$locus = db$gene
db[ !is.na(db$motif) ,]$locus = paste0(db[ !is.na(db$motif) ,]$gene, " (", db[ !is.na(db$motif) ,]$motif, ")")


df[df$inheritance %in% "dominant",]$inheritance = "longer allele"
df[df$inheritance %in% "recessive",]$inheritance = "both alleles"

#both alleles only for recessive RES

df = df[!(df$inheritance %in% "both alleles" & !(df$locus %in% db$locus)),]

df$inheritance = factor(df$inheritance,  levels = rev(unique(df$inheritance)))

xmax =  max(df[df$logCIU < Inf,]$logCIU)*1.1
xmin =  min(df[df$logCIL > -Inf,]$logCIL)*1.1

df$uid = paste0(df$cohort, "-",df$comparison )
df$uid = gsub("AMP-PD-Case-control", "PD", df$uid)
df$uid = gsub("AMP-PD-Other-control", "PD-like", df$uid)
df$uid = gsub("LBD-Other-control", "LBD", df$uid)

df$cohort = factor(df$cohort,  levels = unique(df$cohort))
df$inheritance = factor(df$inheritance,  levels = (unique(df$inheritance)))
df$uid = factor(df$uid,  levels = rev(unique(df$uid)))

annot = df
#annot$lab = paste0(annot$n1, "/", annot$n2, " ; ", annot$n3, "/", annot$n4) 
annot$perc = round(annot$n1/annot$n2*100, digits =2)
annot$lab = paste0( annot$perc, "% (", annot$n1, ")") 

annot$logOR = 11





#df$comparison = factor(df$comparison, ordered = T, levels = unique(df$comparison))


p1 = ggplot(df, aes(x = logOR, y = uid, colour = uid)) + theme_classic() +
  geom_vline(aes(xintercept = 0), linewidth = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = logCIU, xmin = logCIL, colour = uid), 
                 linewidth = .5, height =.5) +
  geom_point(size=2.6, shape=15) + scale_color_manual(values=c( "#27b51d", "#03a9fc", "darkblue")) +
  #scale_size_manual(values=c(2, 3)) +
  #scale_x_continuous(breaks=seq(-20, 20, 4)) +
  theme(    axis.ticks.y = element_blank(),
            panel.grid.minor.x = element_blank(), 
            panel.grid.major.x = element_blank(), 
            panel.grid.major.y = element_line(linetype = 2, linewidth=0.3, colour="grey70"), 
            legend.position = "bottom",
            legend.title=element_blank(),
            axis.text.y = element_blank(), #element_text(face = 'italic', size=8),
            strip.background =element_rect(fill="white"),
            strip.text.x = element_text(size = 12),
            plot.margin = unit(c(1,5,1,1), "lines"),
            panel.spacing.x = unit(4.5, "lines"),
            panel.spacing.y = unit(0.05, "lines"),
            strip.text.y.left = element_text(angle = 0, vjust =1 ,
                                             hjust= 1, face = "italic")) +
  ylab("") +   xlab("Odds ratio (log2)") + 
  geom_text(data = annot, mapping= aes(label = lab), size=3, col= "black", hjust=0) + 
  facet_grid(inheritance + locus ~ threshold , switch = "y",
             #scales="free",
             space = "free")  +
  coord_cartesian(xlim = c(-9, 9), # This focuses the x-axis on the range of interest
                  clip = 'off') +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank()); p1



pdf("Figure_1_Rafehi.pdf", width=9, height = 11)
print(p1)
dev.off()







