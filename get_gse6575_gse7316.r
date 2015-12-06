# Final Project

source("http://bioconductor.org/biocLite.R")

# load all packages
biocLite("simpleaffy")
biocLite("GEOquery")
biocLite("affy")

biocLite("org.Hs.eg.db")

library(simpleaffy)
library(GEOquery)
library(affy)
library(org.Hs.eg.db)

wd <- "~/Documents/Bioinformatics-MedBio/Final_Project/FinalProjectCode/BMBFinalProject" # for SJ

# ANOVA function:
#get differentially expressed genes with ANOVA (function copied from Manny's SeedDevelopment lecture)
doAnova<-function(expvalues, expgroups) {
    anova.results = anova(lm(as.numeric(expvalues)
                             ~ as.factor(expgroups)))$"Pr(>F)"[1]
    anova.results
}


# Getting GSE6575, "Gene expression in blood of children with autism spectrum disorder"

#download.file('http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE6575&format=file', 'gse6575.tar', mode = 'wb')
untar('gse6575.tar', exdir = 'gse6575') #odd error while downloading - check if files are corrupted?

#normalizing the data set
setwd('gse6575')
raw_6575 <- ReadAffy(compress = TRUE, verbose = TRUE)
norm_6575 <- gcrma(raw_6575)
exprs_6575 <- exprs(norm_6575)

#get the experimental groups
geo_6575 <- getGEO("GSE6575", destdir=getwd())
exprs_geo_6575 <- exprs(geo_6575[[1]])
file_designations = colnames(exprs_geo_6575)
design_6575 = vector(length=length(file_designations))
for (i in 1:length(file_designations)) {
  gsm = getGEO(file_designations[i])
  design_6575[i] = Meta(gsm)$characteristics_ch1
}

# Perform ANOVA
anova_6575 <- apply(exprs_6575, 1, doAnova, design_6575)

#########################################

#getting GSE7329, "Gene expression profiles of lymphoblastoid cells (autism)"

# it appears that when you dowload GSE7316_RAW.tar from GEO, you actually get a bunch of txt.gz files, which are already normalized. 
# So these are not 'raw' files. Let's skip getting the data / normalizing and just jump ahead with GEOquery.

#setwd(wd)

geo_7329 <- getGEO('GSE7329') # data might already be in log format, but check with GEO database annotations - SJ
exprs_geo_7329 <- exprs(geo_7329[[1]])


#get the experimental conditions for ANOVA?
file_designations_2 = colnames(exprs_geo_7329)
design_7329 = c(rep('autism_dup15',7),rep('autism_fmr1',8),rep('control',15))
for (i in 1:length(file_designations_2)) {
  gsm = getGEO(file_designations_2[i])
  design_7329[i] = strsplit(Meta(gsm)$title, split=" ")[[1]][1]
}

#get differentially expressed genes with ANOVA (function copied from Manny's SeedDevelopment lecture)
doAnova<-function(expvalues, expgroups) {
  anova.results = anova(lm(as.numeric(expvalues)
                           ~ as.factor(expgroups)))$"Pr(>F)"[1]
  anova.results
}

anova_7329 <- apply(exprs_geo_7329, 1, doAnova, design_7329) #returns ~44,000 elements
anova_7329.fdr = p.adjust(anova_7329, method="BH") # correct for multiple hypothesis testing (FDR)
anova_7329.bfn = p.adjust(anova_7329, method="bonferroni")

test <- which(anova_7329 < .05)
length(test) # ~10,000 elements with signficant p-values, let's see if we can narrow this down somehow
# should p-values be adjusted for multiple-hypothesis testing? - SJ

# check adjusted p-value - SJ
test.fdr <- anova_7329.fdr[anova_7329.fdr < 0.05]
test.bfn <- anova_7329.bfn[anova_7329.bfn < 0.05]
length(test.fdr) # approximately 3882
length(test.bfn) # approximately 149

# get gene names - SJ
genes.fdr <- names(test.fdr)
genes.bfn <- names(test.bfn)
# what are the gene_ids?  Can we convert them to genbank symbols? - SJ
gpl1708 <- read.delim("~/Documents/Bioinformatics-MedBio/Final_Project/FinalProjectCode/BMBFinalProject/GPL1708-20418.txt", row.names=1, comment.char="#")
gpl1708[gpl1708$GENE %in% as.integer(genes.bfn),]

in.gpl <- genes.bfn[as.integer(genes.bfn) %in% gpl1708$GENE]
genes.final <- in.gpl[!is.na(in.gpl)]

# GO term enrichment
final.table <- gpl1708[gpl1708$GENE %in% genes.final, c("GENE","GENE_SYMBOL", "GO_ID")]
p.vals.final <- test.bfn[genes.final]

org.Hs.eg.db



# get log fold change using MKs getRange function - SJ
getRange <- function(expvalues, expgroups) {
    range.results = range(tapply(expvalues, expgroups, mean))
    rangediff = range.results[2]-range.results[1]
}

range.results = apply(exprs_geo_7329, 1, getRange, gse5634_expdesign) # assume exprs_geo_7329 is already log normalized
# finally get the results
diffexp_gse5634 = logexprs_gse5634[range.results > 5 & anova.results.fdr
                                   < 0.05,]


#####
# get gse25507
if(!file.exists("gse225507.tar"))
download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25507&format=file", destfile="gse25507.tar")
untar("gse25507.tar", exdir="gse25507")

setwd("gse25507")

#begin processing data
setwd(paste0(wd,"/gse25507"))
raw_gse25507 <- ReadAffy(compress=TRUE, verbose=TRUE)
norm_gse25507 <- gcrma(gse25507.raw)
exprs_gse25507 <- exprs(gse25507.norm)


geo25507 <- getGEO("GSE25507", destdir=getwd())
exprs_geo25507<- exprs(geo25507[[1]])
file_designations.25507 = colnames(gse25507.exprs)
design_25507 = vector(length=length(file_designations.25507))
for (i in 1:length(file_designations.25507)) {
    filename = strsplit(file_designations.25507[i], split="[.]")[[1]][1]
    print(filename)
    gsm = getGEO(GEO=filename)
    design_25507[i] = paste(Meta(gsm)$characteristics_ch1, collapse=" ")
}

# Perform ANOVA
anova_gse25507 <- apply(exprs_gse25507, 1, doAnova, design_25507)

# Identify significant genes, and also get log fold change


