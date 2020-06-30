

##############################  GTEx data process   ################################

library(data.table)
gtex_data = data.frame(fread("/n/groups/price/kushal/Robocov/data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz"))
gene_names = gtex_data$Description

gtex_mat = gtex_data[,-(1:2)]


library(readxl)
gtex_samples = read.delim2("/n/groups/price/kushal/Robocov/data/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt")

gtex_data2 = gtex_data[match(unique(gene_names), gene_names), ]
gtex_mat2 = gtex_data2[, -(1:2)]

samp_pheno = read.delim2("/n/groups/price/kushal/Robocov/data/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
samp_pheno$NEW_ID = gsub("-", ".", samp_pheno$SAMPID)

tissue_labels = samp_pheno$SMTSD[match(colnames(gtex_mat2), samp_pheno$NEW_ID)]

person_labels = unlist(lapply(strsplit(colnames(gtex_mat2), "[.]"), function(x) return(paste0(x[1:2], collapse = "-"))))

gtex_samples2 = gtex_samples[-980,]


covid19_genes = read.delim("/n/groups/price/kushal/singlecellLDSC/data/Gene_Scores/Modules/healthy/celltype_enriched2/corona/ACE2..TMPRSS2._s.txt",
                           header=F)

idx = match(covid19_genes[,1], gtex_data2$Description)
idx2 = idx[which(!is.na(idx))]
gtex_data3 = gtex_data2[idx2,]
gtex_mat3 = gtex_data3[, -(1:2)]


comp = as.numeric(gtex_mat3[2,])
comp[comp==0] = -999
tmp = xtabs(comp ~ person_labels + tissue_labels)
tmp[tmp == 0] = NA
tmp[tmp == -999] = 0
tmp = as.data.frame.matrix(tmp)
save(tmp, file = "/n/groups/price/kushal/singlecellLDSC/data/ACE2_expression.rda")






