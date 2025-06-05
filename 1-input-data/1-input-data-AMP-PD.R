res = read.delim("raw_input_data/all_results_EH5.txt")

#flag controls that have family history of PD
#leave controls that are carriers for known PD SNPs

res$remove = "keep"
res[(res$cases == "Control" & res$history %in% c("both", "family_history_only")) ,]$remove = "remove_controls"

res$group = paste0(res$cases, "-", res$history)

table(res$history, res$cohort, res$cases)

#flag related
related = read.delim("/vast/projects/bahlo_amppd/str_analysis/flag_related/results_MWVC.txt", header=F)
colnames(related) = "ID"

removed_related = res[res$sample %in% related$ID,]

removed_related_summary = as.data.frame( table(removed_related[!duplicated(removed_related$sample),]$history,
      removed_related[!duplicated(removed_related$sample),]$cohort,
      removed_related[!duplicated(removed_related$sample),]$cases))

write.table(removed_related_summary, file="remove_related_summary.txt", sep="\t", quote = F, row.names = F)

res[(res$sample %in% related$ID),]$remove = "remove_related"



#HAVE NOT REMOVED ANY RELATED PEOPLE BETWEEN COHORTS FOR THIS ANALYSIS. THIS IS ACTUALLY VERY INTERESTING - WHY DO SOME RELATED PEOPLE HAVE PD AND SOME HAVE LBD?
#can add do a specific analysis of the case-case LBD-PD people. for now, leave in the study since we arent collapsing the results 



#Whole cohort
unique = res[!duplicated(res$sample),]
cohort_size = table(unique$cohort, unique$cases)

write.table(cohort_size, file="cohort_size_summary_AMPPD_unfiltered.txt", quote = F, sep = "\t")

cohort_size = read.delim("cohort_size_summary_AMPPD_unfiltered.txt")
cohort_size$total = rowSums(cohort_size)
write.table(cohort_size, file="cohort_size_summary_AMPPD_unfiltered.txt", quote = F, sep = "\t")


#Keep only cohort
unique = unique[unique$remove %in% "keep",]
cohort_size = table(unique$cohort, unique$cases)

write.table(cohort_size, file="cohort_size_summary_AMPPD_filtered.txt", quote = F, sep = "\t")

cohort_size = read.delim("cohort_size_summary_AMPPD_filtered.txt")
cohort_size$total = rowSums(cohort_size)
write.table(cohort_size, file="cohort_size_summary_AMPPD_filtered.txt", quote = F, sep = "\t")



#add in more metadata

metadata = read.delim("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/str_meta_analysis/all_input_data/merged_EH_data_kmer_v2/all_metadata_expanded.txt")
metadata$sex = gsub("Female", "female", metadata$sex)
metadata$sex = gsub("Male", "male", metadata$sex)

res = merge(res, metadata[,c(1,4,5,7)], by="sample", all.x = T)

pca = read.delim("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/ancestry/version2/PC_table_with_estimated_superpop_PC123.txt")
res = merge(res, pca[,c(1, 23)], by="sample", all.x = T)

#switch allele 1 and 2 for males on chrX

db = read.delim("exSTRa_hg38.txt", skip=1)

chrX = db[db$chrom == "chrX",]$locus
head(res)


res[is.na(res$rep2),]$rep2 = 0

res$chrX_switch_allele = ""
res[res$locus %in% chrX & res$rep2 == 0 & res$sex %in% "male",]$chrX_switch_allele = "yes"
res[res$chrX_switch_allele == "yes",]$rep2 = res[res$chrX_switch_allele == "yes",]$rep1

res[res$chrX_switch_allele == "yes",]$rep1 = 0
res$chrX_switch_allele  = NULL





#Remove duplicated individuals:

dups = read.delim("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/metadata/v3/amp_pd_participant_wgs_duplicates.csv", sep=",")
pair1 = intersect(res$sample,dups$participant_id)
pair2 = intersect(res$sample,dups$duplicate_sample_id)

#there are no duplicated samples in the WGS data


res = res[!is.na(res$sample),]

write.table(res, file="all_results_AMPPD_EH5_full.txt", row.names = F, quote = F, sep = "\t")

##Add column to label expansions

res = read.delim("all_results_AMPPD_EH5_full.txt")

db = read.delim("exSTRa_hg38.txt", skip=1)
db$aff_low  = as.numeric(db$aff_low )

loci_keep = read.delim("./loci_keep.txt")
locus_list = loci_keep[loci_keep$keep == "keep",]$locus

res = res[!is.na(res$sample),]

if(any(is.na(res$rep2))){ res[is.na(res$rep2),]$rep2 = 0 }

res$rep1_expanded = "not_expanded"
res$rep2_expanded = "not_expanded"

res = res[res$locus %in% locus_list,]

for(locus in locus_list){
  unstable_low = db[db$locus == locus, ]$aff_low
  
  if(nrow(res[res$rep1 >= unstable_low & res$locus %in% locus,]) > 0){
    res[res$rep1 >= unstable_low & res$locus %in% locus,]$rep1_expanded = 'expanded'
  }
  if(nrow(res[res$rep2 >= unstable_low & res$locus %in% locus,]) > 0){
    res[res$rep2 >= unstable_low & res$locus %in% locus,]$rep2_expanded = 'expanded'
  }

}

write.table(res, file="all_results_AMPPD_EH5_keep_loci_expanded_labelled.txt", sep="\t", quote = F, row.names = F)




##Add column to label expansions - for intermediate threshold alleles

res = read.delim("all_results_AMPPD_EH5_keep_loci_expanded_labelled.txt")

db = read.delim("exSTRa_hg38.txt", skip=1)
db = db[db$has_intermediate == "yes",]
db$incomplete_low  = as.numeric(db$incomplete_low )

loci_keep = read.delim("/stornext/Bioinf/data/lab_bahlo/projects/dementia/PD/AMP_PD/str_meta_analysis/key_input_files/loci_keep.txt")
locus_list = loci_keep[loci_keep$keep == "keep",]$locus
locus_list = intersect(locus_list, db$locus)

res$has_intermediate = "no"
res[res$locus %in% locus_list,]$has_intermediate = "yes"

table( res$locus, res$has_intermediate)
res$rep1_expanded_intermediate = "NA"
res$rep2_expanded_intermediate = "NA"
res[res$has_intermediate == "yes",]$rep1_expanded_intermediate = "not_expanded"
res[res$has_intermediate == "yes",]$rep2_expanded_intermediate = "not_expanded"


for(locus in locus_list){
  inter_low = db[db$locus == locus, ]$incomplete_low
  
  
  if(nrow(res[res$rep1 >= inter_low & res$locus %in% locus,]) > 0){
    res[res$rep1 >= inter_low & res$locus %in% locus,]$rep1_expanded_intermediate = 'expanded'
  }
  if(nrow(res[res$rep2 >= inter_low & res$locus %in% locus,]) > 0){
    res[res$rep2 >= inter_low & res$locus %in% locus,]$rep2_expanded_intermediate = 'expanded'
  }

}

write.table(res, file="all_results_AMPPD_EH5_keep_loci_expanded_labelled.txt", sep="\t", quote = F, row.names = F)






