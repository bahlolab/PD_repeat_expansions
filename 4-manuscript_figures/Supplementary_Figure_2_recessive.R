## run from command line in the directory where the script is located

library(gdata)
library(cowplot)
library(grid)
library(gridExtra)
#library(egg)
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

#Load AMP-PD
res = read.delim("../input_data/all_cohorts_FINAL_EH5_keep_loci_expanded_labelled_annotate_interruptions.txt")

res = res[!(res$locus %in% "FECD3_TCF4"),]
res = res[res$remove %in% "keep",]
res = res[res$ancestry %in% "EUR",]
res = res[!(res$cohort %in% "UKBB-PD"),]
res = res[!(res$cohort %in% "gnomad"),]

#res$main_cohort = "AMP-PD"
#res[res$cohort %in% "LBD",]$main_cohort = "LBD"

colnames(res)
#res = res[-c(5:14, 21, 22 )]


#res[is.na(res$ancestry),]$ancestry = "ND"

db = read.delim("../input_data/exSTRa_hg38.txt", skip=1)

res$cases = gsub("Control", "control", res$cases)
res$cases = gsub("Case", "PD", res$cases)
res[res$cohort %in% "LBD",]$cases = gsub("Other", "LBD", res[res$cohort %in% "LBD",]$cases)
res$cases = gsub("Other", "PD-like", res$cases)

res$cases = factor(res$cases , levels = c("control", "PD", "PD-like", "LBD"))


res$cohort = gsub("AMP-PD", "PD", res$cohort)
res$cohort = gsub("LBD", "LBD", res$cohort)
res$cohort = factor(res$cohort, levels = c("PD", "LBD"))


my_colors = c("#8D8D8D", "darkblue", "#03a9fc", "#27b51d")


db$aff_low  = as.numeric(db$aff_low )
db$norm_up  = as.numeric(db$norm_up )
loci = as.data.frame(unique(res$locus))
colnames(loci) = "locus"
loci = merge(loci, db, by="locus", all.x=T )


max_re = res[,c(3,6)]

max_re1 = max_re %>% 
  group_by(gene) %>%
  summarise(max(rep2))

max_re1$group = "<75"
max_re1[max_re1$`max(rep2)` > 75, ]$group = ">75"
max_re1[max_re1$`max(rep2)` > 1000,]$group = ">1000"

loci_keep = read.delim("../input_data/loci_keep.txt")


all_pops = unique(c(res$ancestry))

c25 <- c(
  "dodgerblue2", "#E31A1C",  "green4",   "#6A3D9A", "#FF7F00",  "gold1",  "skyblue2", "#FB9A99",   "palegreen2",  "#CAB2D6",
  "#FDBF6F",  "gray70", "khaki2",  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",  "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

ancestry_cols = c25[1:13]
names(ancestry_cols) = all_pops


#annotate with gene name
gene_locus_map = res[!duplicated(res$locus),][,c(2,3)]


size_label = 4
size_legend=20
size_main = 20
size_tag=2.5
size_points=5

loci[is.na(loci$incomplete_low),]$incomplete_low = -1000


res$interruptions = "not interrupted"
res[ !is.na(res$interruption_full) ,]$interruptions = "interrupted"
res[ !is.na(res$interruption_int) ,]$interruptions = "interrupted"

## Recessive disorders only
db = read.delim("../input_data/exSTRa_hg38.txt", skip=1)

db = db[db$inheritance %in% c("AR", "XLR"),]
db = db[db$locus %in% loci_keep[loci_keep$keep %in% "keep",]$locus ,]

res1 = res[res$gene %in% c(unique(db$gene), "RFC1_AAGGG", "RFC1_ACAGG"),]

unique(res1$locus)

max_re = res[,c(2,6)]

max_re1 = max_re %>% 
  group_by(locus) %>%
  summarise(max(rep2))
plot_EH = list()

res = res[res$cohort %in% c("PD", "LBD"),]

for (plot_locus in unique(res1$locus)) {
  
  db1 <- loci[loci$locus == plot_locus, ]
  stable_high = db1$norm_up 
  unstable_low = db1$aff_low
  incomplete_low = db1$incomplete_low 
  
  max_re_locus = max_re1[max_re1$locus == plot_locus,]$`max(rep2)`
  upper_lim = max_re_locus*1.4
  
  x = res1[res1$locus == plot_locus,]
  plot_max = max(x$rep2)
  x$index =1
  
  
  
  plot_EH[[plot_locus]] = ggplot(x, 
                                 aes(y = rep2, x = rep1, colour=cases, 
                                     #shape=history, 
                                     alpha=0.6)) + 
    theme_bw() + geom_point(size=size_points) + 
    scale_color_manual(values = my_colors) +
    scale_x_continuous(limits=c(0, upper_lim)) + 
    scale_y_continuous(limits=c(0, upper_lim)) + 
    ylab("longer allele (number of repeats)") + 
    xlab("shorter allele (number of repeats)") + 
    ggtitle( x$gene) + 
    #geom_text(hjust = -0.15, alpha=1, fontface = "bold", show.legend = FALSE) +
    theme(legend.position="bottom", legend.box = "horizontal",
          legend.key.size = unit(0.1, "lines"),
          legend.title=element_blank(), 
          legend.text  = element_text(size = size_legend)
    )       + guides(alpha = "none") +
    theme(text = element_text(size=size_main)) +
    geom_vline(xintercept = unstable_low, col="#D14E4B", linetype=2) +
    geom_hline(yintercept = unstable_low, col="#D14E4B", linetype=2) +
    geom_vline(xintercept = incomplete_low, col="#D19592", linetype=2) +
    geom_hline(yintercept = incomplete_low, col="#D19592", linetype=2) +
    geom_vline(xintercept = stable_high, col="#2528C7",  linetype=2) +
    geom_hline(yintercept = stable_high, col="#2528C7", linetype=2) + 
    facet_grid(. ~ cohort) 
  
  
}


plot_EH


library(ggpubr)

p1 = ggarrange(plot_EH[[1]],plot_EH[[2]], plot_EH[[3]],plot_EH[[4]],plot_EH[[5]],
               ncol = 2, common.legend = T, nrow = 3,
               legend = "bottom")

pdf("Supplementary_Figure_2_recessive.pdf", width=16, height = 16)
print(p1)
dev.off()


