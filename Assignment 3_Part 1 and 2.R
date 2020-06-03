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

#Question 5: Make a histogram plot of the mean values in png format and paste it into your report.
hist(x$Mean1)
hist(x$Mean1, breaks = 20, border="red", col = "green")

#Question 6: Import this csv file into an R object. What are the column names?
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part1_files/growth_data.csv", destfile = "growth_data.csv")
y<- read.csv("growth_data.csv", header = TRUE)
head(y)
str(y)
#The column names are Site, TreeID, Circumf_2004_cm, Circumf_2009_cm, Circumf_2014_cm and Circumf_2019_cm.

#Question 7: Calculate the mean and standard deviation of tree circumference at the start and end of the study at both sites.
mean(y$Circumf_2004_cm)
subset(y, Site=="northeast")
head(y)
str(y)
subset(y, Site=="southwest")
head(y)
str(y)
ne <- subset(y, Site=="northeast")
head(ne)
str(ne)
sw <- subset(y, Site=="southwest")
head(sw)
str(sw)
mean(ne$Circumf_2004_cm)
mean(sw$Circumf_2004_cm)
mean(ne$Circumf_2019_cm)
mean(sw$Circumf_2019_cm)
sd(sw$Circumf_2019_cm)
sd(ne$Circumf_2019_cm)
sd(sw$Circumf_2004_cm)
sd(ne$Circumf_2004_cm)
head(ne)
head(sw)

# Question 8: Make a box plot of tree circumference at the start and end of the study at both sites.
boxplot(ne$Circumf_2004_cm,ne$Circumf_2019_cm)
boxplot(ne$Circumf_2004_cm,ne$Circumf_2019_cm,sw$Circumf_2004_cm,sw$Circumf_2019_cm,names= c("ne2004","ne2019","sw2004","sw2019"), ylab="Circumference (cm)", xlab="site and years", main="Growth at plantation site")






#Assignment 3 Part 2

#downloading library
library("seqinr")
library("R.utils")
library("rBLAST")
library("ape")
library("ORFik")
library("Biostrings")

#Question 1: Download the whole set of E. coli gene DNA sequences and use gunzip to decompress. Use themakeblast() function to create a blast database. How many sequences are present in the E.coli set?
download.file("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-42/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655/cds/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz", destfile = "Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")


R.utils::gunzip("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa.gz")
makeblastdb("Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.cds.all.fa",dbtype="nucl","-parse_seqids")
# 4140 sequences are present in the E.coli set

#Question 2: Download the sample fasta sequences and read them in as above. For your allocated sequence,determine the length (in bp) and the proportion of GC bases.
download.file("https://raw.githubusercontent.com/markziemann/SLE712_files/master/bioinfo_asst3_part2_files/sample.fa", destfile= "sample.fa")
s <- read.fasta("sample.fa")
head(s)
str(s)
#Allocated sequence is "25"
c <- s$`25`
head(c)
str(c)

seqinr::getLength(c)

seqinr::GC(c)
#length for sequence"25" is 1047 and proportion of GC bases is 0.5635148
