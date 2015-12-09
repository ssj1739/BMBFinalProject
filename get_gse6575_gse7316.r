# Final Project
# Sidharth Jain, Warren Ersly, Aisha Dar
source("http://bioconductor.org/biocLite.R")

# load all packages
packages <- c("simpleaffy", 
              "GEOquery", 
              "affy",
              "org.Hs.eg.db",
              "GOstats",
              "GO.db",
              "hgu133plus2.db",
              "cluster",
              "gplots",
              "hgug4112a.db",
              "org.Hs.eg.db")

for(pkg in packages){
    if(!pkg %in% installed.packages()[,1])
        biocLite(pkg)
    require(pkg, character.only=T)
}

# set project directory
wd <- "~/Documents/Bioinformatics-MedBio/Final_Project/FinalProjectCode/BMBFinalProject" # for SJ
# go to project directory (wd)
setwd(wd)

### FUNCTIONS ###

# initialize ANOVA function:
# get differentially expressed genes with ANOVA (function copied from Manny's SeedDevelopment lecture)
doAnova<-function(expvalues, expgroups) {
    anova.results = anova(lm(as.numeric(expvalues)
                             ~ as.factor(expgroups)))$"Pr(>F)"[1]
    return(anova.results)
}

# initialize getRange function (taken from Manny's SeedDevelopment lecture)
getRange <- function(expvalues, expgroups) {
    range.results = range(tapply(expvalues, expgroups, mean))
    rangediff = range.results[2]-range.results[1]
}

#######################################
### BEGIN ANALYSIS ###
# GSE6575 #
setwd(wd)
# Getting GSE6575, "Gene expression in blood of children with autism spectrum disorder"
if(!dir.exists('gse6575')){
    download.file('http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE6575&format=file', 'gse6575.tar', mode = 'wb')
    untar('gse6575.tar', exdir = 'gse6575') #odd error while downloading - check if files are corrupted?
}
setwd('gse6575')

#get the experimental groups
geo_6575 <- getGEO("GSE6575", destdir=getwd())
exprs_geo_6575 <- exprs(geo_6575[[1]])
file_designations_6575 = colnames(exprs_geo_6575)
design_6575 = character(length=length(file_designations))
for (i in 1:length(file_designations)) {
  gsm = getGEO(file_designations[i])
  design_6575[i] = Meta(gsm)$characteristics_ch1
}

# load and normalize data set
covdesc_6575 <- cbind(file_designations_6575, design_6575)
colnames(covdesc_6575) <- c("", "treatment")
write.table(x=covdesc_6575, file="covdesc_6575", sep="   ", row.names=F, quote=F)


files_6575 <- list.files(pattern="CEL.gz")
raw_6575 <- read.affy("covdesc_6575", verbose = TRUE)
norm_6575 <- mas5(raw_6575)
exprs_6575 <- exprs(norm_6575)
log_exprs_6575 <- log2(exprs_6575)


# Perform ANOVA
# full design - all groups accounted for in anova
anova_6575 <- apply(log_exprs_6575, 1, doAnova, design_6575)
fdr_6575 <- p.adjust(anova_6575, method="BH")
bfn_6575 <- p.adjust(anova_6575, method="bonferroni")

# determining significance with log fold change and p-value cutoffs/corrections
range_6575 <- apply(log_exprs_6575, 1, getRange, design_6575) # full

# filter for log fold change and p-value cutoff.  ###!!! NOTE: may want to change cutoff, as logfold change cutoff is not met in full design
diffexp_6575 = log_exprs_6575[fdr_6575 < 0.001 & range_6575 > 1,]
sig_genes_6575 <- rownames(diffexp_6575)
length(sig_genes_6575)

# Annotate genes and perform GO enrichment
universe_6575 <- rownames(exprs_6575)
universe_eid_6575 <- unlist(mget(universe_6575, hgu133plus2ENTREZID))
sig_genes_eid_6575 <- unlist(mget(sig_genes_6575, hgu133plus2ENTREZID))
sig_genes_names <- unlist(mget(sig_genes_6575, hgu133plus2GENENAME))

params_6575=new("GOHyperGParams", geneIds=sig_genes_eid_6575,
           universeGeneIds=universe_eid_6575, annotation="hgu133plus2", ontology="BP",
           pvalueCutoff=0.001, conditional=FALSE, testDirection="over")
overRepresented_6575=hyperGTest(params)

summary(overRepresented_6575)

# clustering/create heat plot
expmatrix_6575 <- exprs_gse6575[gene_list_6575[!is.na(gene_list_6575)],]
expmatrix_cor_6575<-cor(t(expmatrix_6575))
expmatrix_dist_6575<-as.dist(1-expmatrix_cor_6575)

# cluster
expmatrix_hclust_6575<-hclust(
    expmatrix_dist_25507, method="ave"
)

expmatrix_hclust_groups_25507 = cutree(
    expmatrix_hclust_25507, k=2
)

names(
    which(expmatrix_hclust_groups_25507==1)
) ->group1_hclust_25507

names(
    which(expmatrix_hclust_groups_25507==2)
) ->group2_hclust_25507
#pdf("boxplot1_25507.pdf")
boxplot(pair_25507@means[group1_hclust_25507,])
#dev.off()
#pdf("boxplot2_25507.pdf")
boxplot(pair_25507@means[group2_hclust_25507,])
#dev.off()
#silhoutte plot
sil_25507 <- silhouette(expmatrix_hclust_groups_25507, dist=expmatrix_dist_25507)
#pdf("silhoutte_25507.pdf")
plot(sil_25507)
#dev.off()


# create a heatmap
library(gplots)
hclust2 <- function(x, method="average", ...) {
    hclust(x, method=method, ...)
}
dist2 <- function(x, ...) {
    as.dist(1-cor(t(x), method="pearson"))
}
#pdf("heatmap_25507.pdf")
heatmap.2(expmatrix_25507,
          col=redgreen(75),
          hclustfun=hclust2,
          distfun=dist2,
          scale="row",
          cexCol=0.6,
          Colv=TRUE,
          sepcolor="black",
          dendrogram="both",
          key=TRUE,
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow = 0.4
)
dev.off()

#########################################
# GSE7329 #

#getting GSE7329, "Gene expression profiles of lymphoblastoid cells (autism)"

# it appears that when you dowload GSE7316_RAW.tar from GEO, you actually get a bunch of txt.gz files, which are already normalized. 
# So these are not 'raw' files. Let's skip getting the data / normalizing and just jump ahead with GEOquery.

setwd(wd)

# found raw files as text
if(!dir.exists('gse7329')){
    download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE7329&format=file", "gse7329.tar", mode = 'wb')
    untar('gse7329.tar', exdir = 'gse7329') #odd error while downloading - check if files are corrupted?
}
setwd("gse7329")
files_7329 <- as.character(list.files(pattern="GSM", full.names=F))

#get the experimental conditions for ANOVA?
file_designations_2 = colnames(exprs_geo_7329)
design_7329 = c(rep('autism_dup15',7),rep('autism_fmr1',8),rep('control',15))
design_A_7329 = character()
for (i in 1:length(file_designations_2)) {
    gsm_7329 = getGEO(file_designations_2[i], destdir=getwd())
    design_A_7329[i] = strsplit(Meta(gsm_7329)$title, split=" ")[[1]][1]
}

# covdesc_7329 <- cbind(files_7329, design_7329)
# colnames(covdesc_7329) <- c("", "disease")
# write.table(x=covdesc_7329, file="covdesc_7329", sep="   ", row.names=F, quote=F)
# 
# Can't use read.affy, not a CEL format file.  File is already normalized.
# read.affy("covdesc_7329", compress=T, verbose=T)

geo_7329 <- getGEO('GSE7329', destdir=getwd()) # data is already log normalized based on cy3-cy5 ratios.
exprs_geo_7329 <- exprs(geo_7329[[1]])

# get log fold change
range_6575 <- apply(exprs_geo_7329, 1, getRange, design_7329)

# perform anova
anova_7329 <- apply(exprs_geo_7329, 1, doAnova, design_7329) #returns ~44,000 elements
fdr_7329 = p.adjust(anova_7329, method="BH") # correct for multiple hypothesis testing (FDR)
bfn_7329 = p.adjust(anova_7329, method="bonferroni")

# check adjusted p-value - SJ
sig_genes_fdr_7329<- anova_7329.fdr[anova_7329.fdr < 0.05]
sig_genes_bfn_7329 <- anova_7329.bfn[anova_7329.bfn < 0.05]
length(sig_genes_fdr_7329) # approximately 3882
length(sig_genes_bfn_7329) # approximately 149

gene_names_7329_bfn <- names(sig_genes_bfn_7329)
gene_names_7329_fdr <- names(sig_genes_fdr_7329)

# get platform info
gpl1708 <- getGEO(filename='gpl1708.soft', destdir=getwd())

# get gene names - SJ
universe_7329 <- rownames(exprs_geo_7329[3:nrow(exprs_geo_7329),])

# get entrezIDs from gene names
    # first get spot IDs from agilent
conversion_table_1708 <- Table(gpl1708)[,c(1,5,13)]
universe_ESID_7329<-conversion_table_1708[conversion_table_1708[,1] %in% universe_7329,3]
names(sig_genes_fdr_7329) <- conversion_table_1708[conversion_table_1708[,1] %in% names(sig_genes_fdr_7329),3]
names(sig_genes_bfn_7329) <- conversion_table_1708[conversion_table_1708[,1] %in% names(sig_genes_bfn_7329),3]

EID_7329_fdr <- names(sig_genes_fdr_7329[!is.na(sig_genes_fdr_7329)])
EID_7329_bfn <- names(sig_genes_bfn_7329[!is.na(sig_genes_bfn_7329)])

EID_7329_fdr <- EID_7329_fdr[!is.na(EID_7329_fdr) & EID_7329_fdr != ""]
EID_7329_bfn <- EID_7329_bfn[!is.na(EID_7329_bfn) & EID_7329_bfn != ""]

universe_eid_7329 <- unlist(mget(as.character(universe_ESID_7329[universe_ESID_7329 != ""]), org.Hs.egENSEMBL2EG))
sig_genes_fdr_eid_7329 <- unlist(mget(EID_7329_fdr, org.Hs.egENSEMBL2EG))
sig_genes_bfn_eid_7329 <- unlist(mget(EID_7329_bfn, org.Hs.egENSEMBL2EG))
sig_genes_names <- unlist(mget(sig_genes_6575, hgu133plus2GENENAME))

# clustering/create heat plot
expmatrix_7329 <- exprs_geo_7329[[!is.na(gene_list_7329)],]
expmatrix_cor_7329 <- cor(t(expmatrix_7329))
expmatrix_dist_25507<-as.dist(1-expmatrix_cor_25507)

# cluster
expmatrix_hclust_25507<-hclust(
    expmatrix_dist_25507, method="ave"
)

expmatrix_hclust_groups_25507 = cutree(
    expmatrix_hclust_25507, k=2
)

names(
    which(expmatrix_hclust_groups_25507==1)
) ->group1_hclust_25507

names(
    which(expmatrix_hclust_groups_25507==2)
) ->group2_hclust_25507
#pdf("boxplot1_25507.pdf")
boxplot(pair_25507@means[group1_hclust_25507,])
#dev.off()
#pdf("boxplot2_25507.pdf")
boxplot(pair_25507@means[group2_hclust_25507,])
#dev.off()
#silhoutte plot
sil_25507 <- silhouette(expmatrix_hclust_groups_25507, dist=expmatrix_dist_25507)
#pdf("silhoutte_25507.pdf")
plot(sil_25507)
#dev.off()


# create a heatmap
library(gplots)
hclust2 <- function(x, method="average", ...) {
    hclust(x, method=method, ...)
}
dist2 <- function(x, ...) {
    as.dist(1-cor(t(x), method="pearson"))
}
#pdf("heatmap_25507.pdf")
heatmap.2(expmatrix_25507,
          col=redgreen(75),
          hclustfun=hclust2,
          distfun=dist2,
          scale="row",
          cexCol=0.6,
          Colv=TRUE,
          sepcolor="black",
          dendrogram="both",
          key=TRUE,
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow = 0.4
)
dev.off()



##################
# get gse25507
# Download file, if not already downloaded
setwd(wd)
if(!file.exists("gse25507.tar")){
    download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25507&format=file", destfile="gse25507.tar")
    untar("gse25507.tar", exdir="gse25507")
}
setwd("gse25507")

# Get annotations from getGEO, get descriptions for each sample
geo25507 <- getGEO("GSE25507", destdir=getwd())
exprs_geo25507<- exprs(geo25507[[1]])
file_designations_25507 = colnames(gse25507.exprs)
design_25507 = vector(length=length(file_designations.25507))
for (i in 1:length(file_designations.25507)){
    filename = strsplit(file_designations.25507[i], split="[.]")[[1]][1]
    print(filename)
    gsm = getGEO(GEO=filename)
    design_25507[i] = strsplit(Meta(gsm)$characteristics_ch1[2], split=" ")[[1]][2]
}

# get file names
files_25507 <- list.files('gse25507', pattern="CEL.gz", full.names=T)

# create covdesc file
covdesc_25507 <- cbind(file_designations_25507, design_25507)
colnames(covdesc_25507) <- c("", "disease")
write.table(x=covdesc_25507, file="covdesc_25507", sep="   ", row.names=F, quote=F)

raw_25507 <- read.affy("covdesc_25507", verbose = TRUE)
norm_gse25507 <- mas5(raw_25507)
exprs_gse25507 <- exprs(norm_gse25507)

covdesc <- cbind(file_designations_25507, design_25507)
setwd("gse25507")

# Perform pairwise comparison
pair_25507 <- pairwise.comparison(norm_gse25507, group="disease", members=c("autism", "control"), raw_25507)

pair_filt_25507<-pairwise.filter(pair_25507, min.present.no=3, present.by.group=T, fc=log2(1.5), tt=0.05)

gene_list_25507<-rownames(pair_filt_25507@means)

# go term enrichment
universe_25507 <- rownames(exprs_gse25507)
universe_25507_eid <- unlist(mget(universe_25507, hgu133plus2ENTREZID))
gene_eid_25507 <- unlist(mget(gene_list_25507[!is.na(gene_list_25507)], hgu133plus2ENTREZID))
gene_names_25507 <- unlist(mget(gene_list_25507[!is.na(gene_list_25507)], hgu133plus2SYMBOL))

params=new("GOHyperGParams", geneIds=gene_eid_25507,
           universeGeneIds=universe_25507_eid, annotation="hgu133plus2", ontology="BP",
           pvalueCutoff=0.001, conditional=FALSE, testDirection="over")

(overRepresented_25507=hyperGTest(params))

summary(overRepresented_25507)

# clustering/create heat plot
expmatrix_25507 <- exprs_gse25507[gene_list_25507[!is.na(gene_list_25507)],]
expmatrix_cor_25507<-cor(t(expmatrix_25507))
expmatrix_dist_25507<-as.dist(1-expmatrix_cor_25507)

# cluster
expmatrix_hclust_25507<-hclust(
    expmatrix_dist_25507, method="ave"
)

expmatrix_hclust_groups_25507 = cutree(
    expmatrix_hclust_25507, k=2
)

names(
    which(expmatrix_hclust_groups_25507==1)
) ->group1_hclust_25507

names(
    which(expmatrix_hclust_groups_25507==2)
) ->group2_hclust_25507
#pdf("boxplot1_25507.pdf")
boxplot(pair_25507@means[group1_hclust_25507,])
#dev.off()
#pdf("boxplot2_25507.pdf")
boxplot(pair_25507@means[group2_hclust_25507,])
#dev.off()
#silhoutte plot
sil_25507 <- silhouette(expmatrix_hclust_groups_25507, dist=expmatrix_dist_25507)
#pdf("silhoutte_25507.pdf")
plot(sil_25507)
#dev.off()


# create a heatmap
library(gplots)
hclust2 <- function(x, method="average", ...) {
    hclust(x, method=method, ...)
}
dist2 <- function(x, ...) {
    as.dist(1-cor(t(x), method="pearson"))
}
#pdf("heatmap_25507.pdf")
heatmap.2(expmatrix_25507,
          col=redgreen(75),
          hclustfun=hclust2,
          distfun=dist2,
          scale="row",
          cexCol=0.6,
          Colv=TRUE,
          sepcolor="black",
          dendrogram="both",
          key=TRUE,
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow = 0.4
)
dev.off()

# # Perform ANOVA
# anova_gse25507 <- apply(exprs_gse25507, 1, doAnova, design_25507)
# fdr_25507 <- p.adjust(anova_gse25507, method="BH") # correct for multiple hypothesis testing (FDR)
# bfn_25507 <- p.adjust(anova_gse25507, method="bonferroni")
# 
# # Identify significant genes, and also get log fold change
# fdr_25507[fdr_25507 < 0.05]
# bfn_25507[bfn_25507 < 0.05]

