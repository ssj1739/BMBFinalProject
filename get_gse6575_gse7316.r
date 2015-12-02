source("http://bioconductor.org/biocLite.R")
biocLite("simpleaffy")
biocLite("GEOquery")
biocLite("affy")
library(simpleaffy)
library(GEOquery)
library(affy)

wd <- "~/BMBFinalProject"

#getting GSE6575, "Gene expression in blood of children with autism spectrum disorder"
#getting GSE

#download.file('http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE6575&format=file', 'gse6575.tar', mode = 'wb')
untar('gse6575.tar', exdir = 'gse6575') #odd error while downloading - check if files are corrupted?

#get the whole normalized dataset
setwd('gse6575')
raw_6575 <- ReadAffy(compress = TRUE, verbose = TRUE)
norm_6575 <- gcrma(raw_6575)
exprs_6575 <- exprs(norm_6575)

#get the experimental groups
geo_6575 <- getGEO("GSE6575")
exprs_geo_6575 <- exprs(geo_6575[[1]])
file_designations = colnames(exprs_geo_6575)
design_6575 = vector(length=length(file_designations))
for (i in 1:length(file_designations)) {
  gsm = getGEO(file_designations[i])
  design_6575[i] = Meta(gsm)$characteristics_ch1
}

####

#getting GSE7329, "Gene expression profiles of lymphoblastoid cells (autism)"

#it appears that when you dowload GSE7316_RAW.tar from GEO, you actually get a bunch of txt.gz files, which are already normalized. So these are not 'raw' files. Let's skip getting the data / normalizing and just jump ahead with GEOquery.

#setwd('~/r/final_proj_affy_data')

geo_7329 <- getGEO('GSE7329')
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
# shoudl p-values be adjusted for multiple-hypothesis testing?

# check adjusted p-value
test.fdr <- anova_7329.fdr[anova_7329.fdr < 0.05]
test.bfn <- anova_7329.bfn[anova_7329.bfn < 0.05]
length(test.fdr) # approximately 3882
length(test.bfn) # approximately 149

# get gene names
genes.fdr <- names(test.fdr)
genes.bfn <- names(test.bfn)
# what are the gene_ids?  Can we convert them to genbank symbols?


