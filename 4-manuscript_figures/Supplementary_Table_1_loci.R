source('/stornext/Bioinf/data/lab_bahlo/users/rafehi.h/scripts/STRs/all_functions/CanonicalMotif.R')

db = read.delim("../input_data/exSTRa_hg38.txt", skip = 1)
keep = read.delim("../input_data/loci_keep.txt")

db = db[db$locus %in% keep[keep$keep %in% "keep",]$locus,]
db = db[!db$locus %in% "FECD3_TCF4",]
db$coords = paste0(db$chrom, ":", db$hg38_start, "-", db$hg38_end)

res = db[,c("locus", "gene", "long_name", "inheritance", "gene_region", "motif", "coords", "aff_low", "incomplete_low")]
res$locus = gsub("_.*", "", res$locus)
res[res$locus %in% "FGF14", ]$locus = "SCA27B"
res[res$locus %in% "THAP11", ]$locus = "SCA51"
res[res$locus %in% "SCA51", ]$long_name = "Spinocerebellar ataxia 51"

res$motif = CanonicalMotif(res$motif)

res = res[,c(2,1,3:9)]
res$gene_region = gsub("intronic", "intron", res$gene_region)
table(res$gene_region)

colnames(res) = c("gene", "disease abbreviation", "disease full name", "inheritance", "gene region",
                  "canonical motif", "hg38 co-ordinates", "pathogenic threshold", "incomplete penetrance threshold")

write.table(res, file="Supplementary_Table_1_loci.txt", sep="\t", quote = F, row.names = F, col.names = T)
