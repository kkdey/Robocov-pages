

############################  Correlation of gene scores   #####################################

df1 = read.delim("/n/groups/price/kushal/Robocov/data/Gene_Scores/Robospan_mean.txt", header=F)
df2 = read.delim("/n/groups/price/kushal/Robocov/data/Gene_Scores/Corspan_mean.txt", header=F)
df3 = read.delim("/n/groups/price/kushal/Robocov/data/Gene_Scores/pRobospan_mean.txt", header=F)
df4 = read.delim("/n/groups/price/kushal/Robocov/data/Gene_Scores/SEG_GTEx_top10.txt", header=F)



union_genes = Reduce(union, list(df1[,1], df2[,1], df3[,1], df4[,1]))

matt = matrix(0, length(union_genes), 4)
matt[match(df1[,1], union_genes), 1 ] = 1
matt[match(df2[,1], union_genes), 2 ] = 1
matt[match(df3[,1], union_genes), 3 ] = 1
matt[match(df4[,1], union_genes), 4 ] = 1


cor(matt)
