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




summary_data = read.delim("../data/data_summary_counts/data_summary_all_cohort.txt")
summary_data = summary_data[!(summary_data$cohort %in% c("UKBB-PD", "gnomad")),]
summary_data$cases = gsub("Control", "control", summary_data$cases)
summary_data$cases = gsub("Case", "PD", summary_data$cases)
summary_data[summary_data$cohort %in% "LBD",]$cases = gsub("Other", "LBD", summary_data[summary_data$cohort %in% "LBD",]$cases)
summary_data$cases = gsub("Other", "PD-like", summary_data$cases)

summary_data = summary_data[summary_data$ancestry %in% "EUR",]
summary_data = summary_data[!(summary_data$cohort %in% "LBD" & summary_data$cases %in% c("PD", "PD-like")),]
summary_data$label = paste0(summary_data$cases, ": ", summary_data$count_DOM, " (", summary_data$percent_DOM, "%)")

summary_data$cohort = gsub("AMP-PD", "PD", summary_data$cohort)
summary_data$cohort = gsub("LBD", "LBD", summary_data$cohort)
summary_data$cohort = factor(summary_data$cohort, levels = c("PD", "LBD"))

#annotate with gene name
gene_locus_map = res[!duplicated(res$locus),][,c(2,3)]

summary_data = merge(summary_data, gene_locus_map, by="locus")

size_label = 4
size_legend=25
size_main = 25
size_tag=2.5
size_points=5

loci[is.na(loci$incomplete_low),]$incomplete_low = -1000


res$interruptions = "not interrupted"
res[ !is.na(res$interruption_full) ,]$interruptions = "interrupted"
res[ !is.na(res$interruption_int) ,]$interruptions = "interrupted"

#### GROUP 1

#for(group in unique(max_re1$group)[1]){
group = unique(max_re1$group)[1]

plot_locus = max_re1[max_re1$group %in% group,]$gene
x = res[res$gene %in% plot_locus,]
plot_max = max(x$rep2)
x$index =1

db1 = loci[loci$gene %in% plot_locus, ]
vline = merge(db1[,c(1,5,15,20)], res, by="gene")
vline$uid = paste0(vline$cohort, vline$gene)
vline=vline[!duplicated(vline$uid),]
vline$index = 1
#vline[vline$aff_low > 75,]$aff_low = NA_character_
#vline[vline$incomplete_low > 75,]$incomplete_low = NA_character_



x1 = summary_data[summary_data$gene %in% plot_locus & summary_data$expanded %in% "yes",]
x1$rep2 =  0
x1$index = 1.6
x2 = x1[x1$cases %in% "control",]
x3 = x1[x1$cases %in% "PD",]
x3$index = 1.3
x4 = x1[x1$cases %in% c("PD-like"),]
x4$index = 1
x5 = x1[x1$cases %in% c("LBD"),]
x5$index = 1.3




plot_swimlane = ggplot(x, aes(x=rep2, y=index, colour=cases)) + 
  geom_jitter(width=0, height = 0.5, alpha=0.5, size=4) + facet_grid(gene ~ cohort) +
  xlab("longer allele") + theme_bw() + 
  scale_color_manual(values = my_colors) +
  ylim(c(0.3,1.7)) + 
  xlim(c(0,75)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.key.size = unit(0.1, "lines"),
        legend.title=element_blank(), 
        legend.text  = element_text(size = size_legend),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect(fill="white") ,
        text = element_text(size=size_main),
        strip.text.y.right = element_text(angle = 0, hjust= 0, face = "italic")) +
  scale_alpha(guide="none") + 
  geom_vline(data=vline, aes(xintercept = aff_low), col="#912323", linetype=2, linewidth=1) +
  geom_vline(data=vline, aes(xintercept = incomplete_low), col="#dba042", linetype=2, linewidth=1) 



#png(paste0("plots/","SUPFIG1_swimlane_",group,".png"), width=3000, height = 4000, res= 300)
#print(plot_swimlane)
#dev.off()



#### GROUP 2

#for(group in unique(max_re1$group)[1]){
group = unique(max_re1$group)[2]

plot_locus = max_re1[max_re1$group %in% group,]$gene
x = res[res$gene %in% plot_locus,]
plot_max = max(x$rep2)
x$index =1

db1 = loci[loci$gene %in% plot_locus, ]
vline = merge(db1[,c(1,5,15,20)], res, by="gene")
vline$uid = paste0(vline$cohort, vline$gene)
vline=vline[!duplicated(vline$uid),]
vline$index = 1
#vline[vline$aff_low > 75,]$aff_low = NA_character_
#vline[vline$incomplete_low > 75,]$incomplete_low = NA_character_

max = max(x$rep2)


plot_swimlane2 = ggplot(x, aes(x=rep2, y=index, colour=cases)) + 
  geom_jitter(width=0, height = 0.5, alpha=0.5, size=4) + facet_grid(gene ~ cohort) +
  xlab("longer allele") + theme_bw() + 
  scale_color_manual(values = my_colors) +
  ylim(c(0.3,1.7)) + 
  xlim(c(0,max*1.1)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.key.size = unit(0.1, "lines"),
        legend.title=element_blank(), 
        legend.text  = element_text(size = size_legend),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect(fill="white") ,
        text = element_text(size=size_main),
        strip.text.y.right = element_text(angle = 0, hjust= 0, face = "italic")) +
  scale_alpha(guide="none") + 
  geom_vline(data=vline, aes(xintercept = aff_low), col="#912323", linetype=2, linewidth=1) +
  geom_vline(data=vline, aes(xintercept = incomplete_low), col="#dba042", linetype=2, linewidth=1) 



#png(paste0("plots/","SUPFIG1_swimlane_",group,".png"), width=2000, height = 4000, res= 300)
#print(plot_swimlane2)
#dev.off()








#### GROUP 3

#for(group in unique(max_re1$group)[1]){
group = unique(max_re1$group)[3]

plot_locus = max_re1[max_re1$group %in% group,]$gene
x = res[res$gene %in% plot_locus,]
plot_max = max(x$rep2)
x$index =1

loci$gene = toupper(loci$gene)

db1 = loci[loci$gene %in% plot_locus, ]
vline = merge(db1[,c(1,5,15,20)], res, by="gene")
vline$uid = paste0(vline$cohort, vline$gene)
vline=vline[!duplicated(vline$uid),]
vline$index = 1
#vline[vline$aff_low > 75,]$aff_low = NA_character_
#vline[vline$incomplete_low > 75,]$incomplete_low = NA_character_

max = max(x$rep2)


plot_swimlane3 = ggplot(x, aes(x=rep2, y=index, colour=cases)) + 
  geom_jitter(width=0, height = 0.5, alpha=0.5, size=4) + facet_grid(gene ~ cohort) +
  xlab("longer allele") + theme_bw() + 
  scale_color_manual(values = my_colors) +
  ylim(c(0.3,1.7)) + 
  xlim(c(0,max*1.1)) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
  theme(legend.position="bottom", legend.box = "horizontal",
        legend.key.size = unit(0.1, "lines"),
        legend.title=element_blank(), 
        legend.text  = element_text(size = size_legend),
        panel.spacing = unit(0, "lines"),
        strip.background = element_rect(fill="white") ,
        text = element_text(size=size_main),
        strip.text.y.right = element_text(angle = 0, hjust= 0, face = "italic")) +
  scale_alpha(guide="none") + 
  geom_vline(data=vline, aes(xintercept = aff_low), col="#912323", linetype=2, linewidth=1) +
  geom_vline(data=vline, aes(xintercept = incomplete_low), col="#dba042", linetype=2, linewidth=1) 



#png(paste0("plots/","SUPFIG1_swimlane_",group,".png"), width=3000, height = 1000, res= 300)
#print(plot_swimlane3)
#dev.off

library(ggpubr)

p1 = ggarrange(plot_swimlane, plot_swimlane2, plot_swimlane3, ncol = 1, common.legend = T,
               legend = "bottom", heights = c(9,5,1.4))

png("Supplementary_Figure_1_swimlane.png", width=3500, height = 6000, res= 300)
print(p1)
dev.off()


