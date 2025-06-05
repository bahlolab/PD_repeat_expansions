## run from command line in the directory where the script is located
set.seed(2)

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

res = res[res$remove %in% "keep",]
res = res[res$ancestry %in% "EUR",]
res = res[!(res$cohort %in% "UKBB-PD"),]

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
res$cohort = gsub("gnomad", "gnomAD", res$cohort)
res$cohort = factor(res$cohort, levels = c("PD", "LBD", "gnomAD"))




my_colors = c("#8D8D8D", "darkblue", "#03a9fc", "#27b51d")


db$aff_low  = as.numeric(db$aff_low )
db$norm_up  = as.numeric(db$norm_up )
loci = as.data.frame(unique(res$locus))
colnames(loci) = "locus"
loci = merge(loci, db, by="locus", all.x=T )


max_re = res[,c(2,6)]

max_re1 = max_re %>% 
  group_by(locus) %>%
  summarise(max(rep2))


loci_keep = read.delim("../input_data/loci_keep.txt")

plot_locus_list = loci_keep[loci_keep$keep == "keep",]$locus
no_hits_after_kmers = setdiff(plot_locus_list, max_re1$locus)
plot_locus_list = intersect(plot_locus_list, max_re1$locus)

all_pops = unique(c(res$ancestry))

c25 <- c(
  "dodgerblue2", "#E31A1C",  "green4",   "#6A3D9A", "#FF7F00",  "gold1",  "skyblue2", "#FB9A99",   "palegreen2",  "#CAB2D6",
  "#FDBF6F",  "gray70", "khaki2",  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",  "darkturquoise", "green1", "yellow4", "yellow3","darkorange4", "brown")

ancestry_cols = c25[1:13]
names(ancestry_cols) = all_pops



#system("mkdir -p other_swimlane/")


summary_data = read.delim("../data/data_summary_counts/data_summary_all_cohort.txt")
summary_data = summary_data[!(summary_data$cohort %in% "UKBB-PD"),]
summary_data$cases = gsub("Control", "control", summary_data$cases)
summary_data$cases = gsub("Case", "PD", summary_data$cases)
summary_data[summary_data$cohort %in% "LBD",]$cases = gsub("Other", "LBD", summary_data[summary_data$cohort %in% "LBD",]$cases)
summary_data$cases = gsub("Other", "PD-like", summary_data$cases)

summary_data = summary_data[summary_data$ancestry %in% "EUR",]
summary_data = summary_data[!(summary_data$cohort %in% "LBD" & summary_data$cases %in% c("PD", "PD-like")),]
summary_data = summary_data[!(summary_data$cohort %in% "gnomad" & summary_data$cases %in% c("PD-like", "LBD", "PD")),]
summary_data$label = paste0(summary_data$cases, ": ", summary_data$count_DOM, " (", summary_data$percent_DOM, "%)")

summary_data$cohort = gsub("AMP-PD", "PD", summary_data$cohort)
summary_data$cohort = gsub("LBD", "LBD", summary_data$cohort)
summary_data$cohort = gsub("gnomad", "gnomAD", summary_data$cohort)
summary_data$cohort = factor(summary_data$cohort, levels = c("PD", "LBD", "gnomAD"))


size_label = 4
size_legend=15
size_main = 16
size_tag=2.5
size_points=5

loci[is.na(loci$incomplete_low),]$incomplete_low = -1000


res$interruptions = "not interrupted"
res[ !is.na(res$interruption_full) ,]$interruptions = "interrupted"
res[ !is.na(res$interruption_int) ,]$interruptions = "interrupted"


plot_EH = list()
for (plot_locus in plot_locus_list[29]) {
  
  db1 <- loci[loci$locus == plot_locus, ]
  stable_high = db1$norm_up 
  unstable_low = db1$aff_low
  incomplete_low = db1$incomplete_low 
  
  max_re_locus = max_re1[max_re1$locus == plot_locus,]$`max(rep2)`
  upper_lim = max_re_locus*1.4
  
  x = res[res$locus == plot_locus,]
  plot_max = max(x$rep2)
  x$index =1
  
  x[x$cohort %in% "gnomAD" & x$rep2 == 34,]$interruptions = "not interrupted"
  x[x$cohort %in% "gnomAD" & x$rep2 == 33,]$interruptions = "interrupted"
  x[x$cohort %in% "gnomAD" & x$rep2 == 33,]$interruption_int = 1
  
  
  
  x1 = summary_data[summary_data$locus %in% plot_locus & summary_data$expanded %in% "yes",]
  x1$rep2 =  8
  x1$index = 1.5
  x2 = x1[x1$cases %in% "control",]
  x3 = x1[x1$cases %in% "PD",]
  x3$index = 1.1
  x4 = x1[x1$cases %in% c("PD-like"),]
  x4$index = 0.7
  x5 = x1[x1$cases %in% c("LBD"),]
  x5$index = 1.1
  
  
  plot_swimlane = ggplot(x, aes(x=rep2, y=index, colour=cases)) + 
    geom_jitter(aes(shape=interruptions), 
                width=0, height = 0.5, alpha=0.5, size=4) + facet_grid(cohort ~ ., switch ="y")  +
    xlab("longer allele") + theme_bw() + 
    scale_color_manual(values = my_colors) +
    scale_shape_manual(values=c(17, 19))+
    xlim(c(8,plot_max)) +  ylim(c(0.3,1.7)) + 
    theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    theme(legend.position="bottom", legend.box = "horizontal",
          legend.key.size = unit(0.1, "lines"),
          legend.title=element_blank(), 
          legend.text  = element_text(size = size_legend),
          panel.spacing = unit(0, "lines"),
          strip.background = element_rect(fill="white") ,
          strip.text.y.right = element_text(angle = 0, hjust= 0),
          text = element_text(size=size_main)) + scale_alpha(guide="none") + 
    ggtitle( x$gene[1]) + 
    scale_x_continuous(breaks = seq(0, 40, by = 2)) + 
    geom_vline(xintercept = unstable_low,col="#912323", linetype=2) +
    geom_vline(xintercept = incomplete_low, col="#dba042", linetype=2) +
    #geom_vline(xintercept = stable_high,  col="#2528C7", linetype=2) + 
    geom_text(data = x2, mapping= aes(label = label), size=size_main/4, col= "black", hjust=0.2) +
    geom_text(data = x3, mapping= aes(label = label), size=size_main/4, col= "black", hjust=0.2) +
    geom_text(data = x4, mapping= aes(label = label), size=size_main/4, col= "black", hjust=0.2) +
    geom_text(data = x5, mapping= aes(label = label), size=size_main/4, col= "black", hjust=0.2) ; plot_swimlane
  
  
  #panel1a = plot_grid(xdens, plot_swimlane, ncol = 1, rel_heights = c(1,3))
 # pdf(paste0("other_swimlane/swimlane_", plot_locus, "_all_cohorts_short_list_panel_with_ancestry.pdf"), height=4, width=8)
  #print(plot_swimlane)
  #dev.off()
  

  
  #plot_grid(plot_EH,  panel1a,  ncol = 2, nrow=1, rel_widths = c(3,1))
  
}






#interuptions
library(ggplot2)
atxn2a = read.delim("../explore_hits/ATXN2/hits_AMP-PD_reviewer_check_SCA2_with_notes-annotated.txt")
atxn2b = read.delim("../explore_hits/ATXN2/hits_AMP-PD_reviewer_check_SCA2_with_notes_top_up-annotated.txt")
atxn2 = rbind(atxn2a, atxn2b)
atxn2 = atxn2[atxn2$ancestry %in% "EUR",]

res = res[res$cohort %in% c("PD", "LBD"),]
all_atxn2 = res[res$gene %in% "ATXN2" ,]
all_atxn2 = as.data.frame(table(all_atxn2$rep2), stringsAsFactors = F)
df = as.data.frame(c(22:39))
colnames(df) = "Var1"

all_atxn2 = merge(all_atxn2, df, by="Var1", all = T)
all_atxn2[is.na(all_atxn2$Freq),]$Freq = 0


#all_atxn2 = all_atxn2[all_atxn2$Freq > 0,]

expanded_atxn2 = res[res$gene %in% "ATXN2" & res$rep2 > 34,]

#atxn2[atxn2$interruption %in% "1,2",]$interruption = 2  

atxn2$interruption = as.numeric(atxn2$interruption)

hist(atxn2$interruption)

x = as.data.frame(table(atxn2$rep2, atxn2$interruption))
y = as.data.frame(table(expanded_atxn2$rep2, expanded_atxn2$interruption_full))

x = rbind(x,y)

x$expanded = "1"
x$Var1 = as.numeric(as.character(x$Var1))
x[x$Var1 >= 34,]$expanded = "2"

library(ggbreak) 

all1 = ggplot(data = all_atxn2, aes(x=Var1, y=Freq)) + geom_col(fill="darkblue") + theme_bw() +
  scale_x_discrete(breaks = seq(19, 39, by = 1)) + scale_y_break(c(250, 500), scales=0.2) +
  xlab("number of repeats") + ylab("number of alleles") + #scale_y_continuous(trans='log2') +
  ggtitle("ATXN2") ; all1

all2 = ggplot(data = all_atxn2, aes(x=Var1, y=Freq)) + geom_col(fill="darkblue") + theme_bw() +
  scale_x_discrete(breaks = seq(19, 39, by = 1)) + #scale_y_break(c(250, 1000), scales=0.2) +
  #scale_y_continuous(breaks = seq(0, 16, by = 2)) + 
  xlab("number of repeats") + ylab("number of alleles") + scale_y_continuous(trans='log2') +
  geom_vline(xintercept = unstable_low+0.5, col="#912323", linetype=2) +
  geom_vline(xintercept = incomplete_low+0.5, col="#dba042", linetype=2) +
  ggtitle("ATXN2 longer allele distribution") ; all2

all_atxn2$Var1 = as.numeric(all_atxn2$Var1)
all3 = ggplot(data = all_atxn2, aes(x=Var1, y=Freq)) + geom_col(fill="darkblue") + theme_bw() +
  scale_x_continuous(breaks = seq(19, 39, by = 1)) + 
  xlab("number of repeats") + ylab("number of alleles") +
  geom_vline(xintercept = unstable_low+0.5, col="#912323", linetype=2) +
  geom_vline(xintercept = incomplete_low+0.5, col="#dba042", linetype=2) +
  theme(text = element_text(size=size_main)) + 
  ggtitle("ATXN2 longer allele distribution") ; all3

#x$Var1 = as.character(x$Var1)
int = ggplot(data = x, aes(x=Var1, y=Freq, fill=Var2)) + geom_col() + theme_bw() +
  xlab("number of repeats") + ylab("number of alleles") + 
  #xlim(21, 39) + 
  scale_x_continuous(breaks = seq(21, 39, by = 1)) + 
  scale_y_continuous(breaks = seq(0, 10, by = 2)) + 
  geom_vline(xintercept = unstable_low+0.5, col="#912323", linetype=2) +
  geom_vline(xintercept = incomplete_low+0.5, col="#dba042", linetype=2) +
  guides(fill=guide_legend(title="number of interruptions", nrow=1),) + 
  scale_fill_brewer(palette = "Greens") + 
  theme(legend.position = "bottom",
        legend.text  = element_text(size = size_legend),
        text = element_text(size=size_main))  + 
  ggtitle("C>T interruptions in the ATXN2 longer alelle") ; int

panel = plot_grid( plot_swimlane, all3, int, ncol = 1,
                   rel_heights = c(1,0.8,1 ), labels = c("A", "B", "C"))

pdf("Figure_2_Rafehi.pdf", height=12, width=8)
print(panel)
dev.off()

png("Figure_2_Rafehi.png", height=2500, width=2500, res=300)
print(panel)
dev.off()



1/700*100

3/569*100
