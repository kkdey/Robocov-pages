

####################################  SEG blood genes   ##############################################

seg_gtex = read.table("/Users/kushaldey/Documents/Robocov-pages/Gene_Sets/SEG_GTEx_top10.txt")
dim(seg_gtex)

gene_names_gtex = as.character(read.table("/Users/kushaldey/Documents/Robocov-pages/data/gene_names_GTEx.txt")[,1])

df = read.table("/Users/kushaldey/Documents/Robocov-pages/output/gene_anno_unique_datefix.txt", header=T)
df2 = df[which(df$id %in% gene_names_gtex), ]


filename = "NMF_run_5_standard"
ll = list.files(paste0("/Users/kushaldey/Documents/Robocov-pages/output/", filename), pattern = ".txt")

for(num in 1:length(ll)){
  test1 = as.character(read.delim(paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/",
                                         filename, "/",
                                         ll[num]),
                                  header = F)[,1])
  cc = as.character(seg_gtex[,1])
  tt = length(intersect(cc, test1))/length(intersect(cc,df2$symbol))/ (length(test1)/length(df2$symbol))
  print(tt)
}


####################################  housekeeping genes   ##############################################

housekeep = read.table("/Users/kushaldey/Documents/Robocov-pages/data/housekeeping_genes.txt")
head(housekeep)

filename = "NMF_run_5_probocov"
ll = list.files(paste0("/Users/kushaldey/Documents/Robocov-pages/output/", filename), pattern = ".txt")

for(num in 1:length(ll)){
  test1 = as.character(read.delim(paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/",
                                         filename, "/",
                                         ll[num]),
                                  header = F)[,1])
  cc = as.character(housekeep[,1])
  tt = length(intersect(cc, test1))/length(intersect(cc,df2$symbol))/ (length(test1)/length(df2$symbol))
  print(tt)
}

####################################  drug targets   ##############################################


drug_targets = data.frame(readxl::read_excel("/Users/kushaldey/Documents/Mouse_Humans/data/gold_standard_drugs_PI.xlsx"))
drug_genes = drug_targets$Gene

filename = "NMF_run_5_standard"
ll = list.files(paste0("/Users/kushaldey/Documents/Robocov-pages/output/", filename), pattern = ".txt")

for(num in 1:length(ll)){
  test1 = as.character(read.delim(paste0("/Users/kushaldey/Documents/Robocov-pages/data/Gene_Scores/",
                                         filename, "/",
                                         ll[num]),
                                  header = F)[,1])
  cc = as.character(drug_genes)
  tt = length(intersect(cc, test1))/length(intersect(cc,df2$symbol))/ (length(test1)/length(df2$symbol))
  print(tt)
}



