#Assignment 3 Part 1
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/gene_expression.tsv", destfile = "gene.expression.tsv")

#Question 1:Read in the file, making the gene accession numbers the row names. Show a table of values for the first six genes. 
x<- read.table("gene.expression.tsv")
x <-read.table("gene.expression.tsv",header= TRUE, row.names = 1)
x[1:7, ]
head(x)
str(x)

#Question 2:Make a new column which is the mean of the other columns. Show a table of values for the first six genes. 
x$Mean1 <- rowMeans(x)
head(x)

#Question 3:  List the 10 genes with the highest mean expression 
order(-x$Mean1)
x[order(-x$Mean1), ]
x_sorted <- x[order(-x$Mean1), ]
head(x_sorted,10)

#Question 4:  Determine the number of genes with a mean <10 
a <- subset(x, Mean1<10) 
head(a)
str(a)

