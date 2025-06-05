system("mkdir alt_thresholds")

#load gnomad data 
gnomad = read.delim("/stornext/Bioinf/data/lab_bahlo/public_datasets/gnomAD3_STRs/gnomAD_STR_genotypes__2022_01_20.tsv")
map = read.delim("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/str_meta_analysis/key_input_files/ExpansionHunter_hg38_gnomad_map.txt")

gnomad$uid = paste(gnomad$Id, gnomad$Motif, sep="_")
#fix alleles so Allele2 used off target reads for FXS_FMR1 and C9orf72
gnomad[gnomad$uid == "C9ORF72_GGCCCC",]$Allele2 = gnomad[gnomad$uid == "C9ORF72_GGCCCC",]$Allele2UsingOfftargetRegions
gnomad[gnomad$uid == "FMR1_CGG",]$Allele2 = gnomad[gnomad$uid == "FMR1_CGG",]$Allele2UsingOfftargetRegions




gnomad$motif1 = gnomad$Motif
gnomad$motif2 = gnomad$Motif
gnomad[grepl("/", gnomad$Motif),]$motif1 = gsub("/.*", "", gnomad[grepl("/", gnomad$Motif),]$Motif)
gnomad[grepl("/", gnomad$Motif),]$motif2 = gsub(".*/", "", gnomad[grepl("/", gnomad$Motif),]$Motif)

gnomad_single = gnomad[!grepl("/", gnomad$Motif),]
gnomad_split1 = gnomad[grepl("/", gnomad$Motif),]
gnomad_split2 = gnomad[grepl("/", gnomad$Motif),]

gnomad_split1$Motif = gnomad_split1$motif1
gnomad_split1$Allele2 = gnomad_split1$Allele1
gnomad_split1$Allele1 = 0

gnomad_split2$Motif   = gnomad_split2$motif2
gnomad_split2$Allele1 = 0

gnomad = rbind(gnomad_single, gnomad_split1, gnomad_split2)
gnomad$motif1 = gnomad$motif2 = NULL




gnomad_cohort_summary = as.data.frame( table(gnomad$Id) )
gnomad_cohort_size = max(gnomad_cohort_summary$Freq)

#rename new UID with new motif
gnomad$uid = paste(gnomad$Id, gnomad$Motif, sep="_")


db = read.delim("exSTRa_hg38.txt", skip=1)
db$aff_low  = as.numeric(db$aff_low )

loci_keep = read.delim("loci_keep.txt")
locus_list = loci_keep[loci_keep$keep == "keep",]$locus

gnomad$rep1_expanded = "not_expanded"
gnomad$rep2_expanded = "not_expanded"

#Map IDs across

gnomad1 = merge(gnomad, map[,c(1,4)], by.x="uid", by.y="gnomad_id")
gnomad = gnomad1

gnomad = gnomad[gnomad$locus %in% locus_list,]

#swap chrX in males from allele 1 to allele 2:

gnomad[gnomad$Chrom == "chrX" & gnomad$Sex == "XY",]$Allele2 = gnomad[gnomad$Chrom == "chrX" & gnomad$Sex == "XY",]$Allele1
gnomad[gnomad$Chrom == "chrX" & gnomad$Sex == "XY",]$Allele1 = 0




for(locus in locus_list){
  unstable_low = db[db$locus == locus, ]$aff_low
  
  if(nrow(gnomad[gnomad$Allele1 >= unstable_low &gnomad$locus %in% locus,]) > 0){
    gnomad[gnomad$Allele1 >= unstable_low &gnomad$locus %in% locus,]$rep1_expanded = 'expanded'
  } 
  
  if(nrow(gnomad[gnomad$Allele2 >= unstable_low &gnomad$locus %in% locus,]) > 0){
    gnomad[gnomad$Allele2 >= unstable_low &gnomad$locus %in% locus,]$rep2_expanded = 'expanded'
  } 
  
  }
  
  
# intermediate
db = db[db$has_intermediate %in% "yes",]

gnomad$has_intermediate = "no"
gnomad[gnomad$locus %in% db$locus,]$has_intermediate = "yes"
gnomad$rep1_expanded_intermediate = "NA"
gnomad$rep2_expanded_intermediate = "NA"
gnomad[gnomad$has_intermediate == "yes",]$rep1_expanded_intermediate = "not_expanded"
gnomad[gnomad$has_intermediate == "yes",]$rep2_expanded_intermediate = "not_expanded"

for(locus in locus_list){
  unstable_low = db[db$locus == locus, ]$incomplete_low
  
  if(nrow(gnomad[gnomad$Allele1 >= unstable_low &gnomad$locus %in% locus,]) > 0){
    gnomad[gnomad$Allele1 >= unstable_low &gnomad$locus %in% locus,]$rep1_expanded_intermediate = 'expanded'
  }  
  
  if(nrow(gnomad[gnomad$Allele2 >= unstable_low &gnomad$locus %in% locus,]) > 0){
    gnomad[gnomad$Allele2 >= unstable_low &gnomad$locus %in% locus,]$rep2_expanded_intermediate = 'expanded'
  }  
  
  }



gnomad_summary = as.data.frame(table(gnomad$locus, gnomad$rep1_expanded))

gnomad$uid = NULL
gnomad = gnomad[,c(24, 1:23, 25:27)]
colnames(gnomad)[2] = "gene"
gnomad$LocusId = NULL
gnomad$cohort = "gnomad"

gnomad = gnomad[,c(21, 1, 14, 15, 27, 10, 9, 11, 22:26)]

colnames(gnomad)[1:8] =c("sample", "locus", "rep1", "rep2", "cohort", "sex","population", "age")

#relabel sex
gnomad$sex = gsub("XY", "male", gnomad$sex)
gnomad$sex = gsub("XX", "female", gnomad$sex)
table(gnomad$sex)


gnomad[gnomad$sample %in% "",]$sample = 1:length(gnomad[gnomad$sample %in% "",]$sample)

gnomad[gnomad$sample %in% "1",]

write.table(gnomad, file="all_results_gnomad_EH5_keep_loci_expanded_labelled.txt", sep="\t", row.names = F, quote = F)

