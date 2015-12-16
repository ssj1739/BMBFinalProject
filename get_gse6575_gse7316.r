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
              "org.Hs.eg.db",
              "Category")

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
design_6575 = character(length=length(file_designations_6575))
for (i in 1:length(file_designations_6575)) {
  gsm = getGEO(file_designations_6575[i])
  design_6575[i] = strsplit(Meta(gsm)$characteristics_ch1, split=" ")[[1]][1]
}
design_full_6575 <- cbind(sapply(file_designations_6575, paste0, ".CEL.gz"), design_6575)

design_6575_pairwise <- design_full_6575[design_full_6575[,2] != "mental",]



# load and normalize data set
if(!file.exists("covdesc")){
    covdesc_6575 <- cbind(design_6575_pairwise[,1], design_6575_pairwise[,2])
    colnames(covdesc_6575) <- c("", "treatment")
    covdesc_final_6575 <- covdesc_6575[covdesc_6575[,2] == "Autism" | covdesc_6575[,2]=="General",] # neglect all mental retardation samples
    write.table(x=covdesc_final_6575, file="covdesc_6575", sep="   ", row.names=F, quote=F)
}
raw_6575 <- read.affy("covdesc_6575", verbose = TRUE)
norm_6575 <- mas5(raw_6575)
exprs_6575 <- exprs(norm_6575)
log_exprs_6575 <- log2(exprs_6575)

# perform pairwise test between autism and control (general pop.)
pair_6575 <- pairwise.comparison(norm_6575, group="treatment", members=c("Autism", "General"), raw_6575)
pair_filt_6575 <- pairwise.filter(pair_6575, min.present.no=3, present.by.group=T, fc=0, tt=0.05)

gene_list_6575 <- rownames(pair_filt_6575@means)
gene_list_6575_narm <- gene_list_6575[!is.na(gene_list_6575)]
length(gene_list_6575_narm) 

# determining significance with log fold change and p-value cutoffs/corrections
range_6575 <- apply(log_exprs_6575, 1, getRange, design_6575_pairwise) # pairwise

# filter for log fold change and p-value cutoff.  ###!!! NOTE: may want to change cutoff, as logfold change cutoff is not met in full design
names_genes_pairfilt_6575 <- names(pair_filt_6575@tt[pair_filt_6575@tt < 0.05])
diffexp_6575 = log_exprs_6575[names_genes_pairfilt_6575[!is.na(names_genes_pairfilt_6575)],]
sig_genes_6575 <- rownames(diffexp_6575)
length(sig_genes_6575)

# Annotate genes and perform GO enrichment
universe_6575 <- rownames(exprs_6575)
universe_eid_6575 <- unlist(mget(universe_6575, hgu133plus2ENTREZID, ifnotfound = NA))
sig_genes_eid_6575 <- unlist(mget(sig_genes_6575, hgu133plus2ENTREZID))
sig_genes_names <- unlist(mget(sig_genes_6575, hgu133plus2GENENAME))

params_6575=new("GOHyperGParams", geneIds=sig_genes_eid_6575,
           universeGeneIds=universe_eid_6575, annotation="hgu133plus2", ontology="BP",
           pvalueCutoff=0.001, conditional=FALSE, testDirection="over")
overRepresented_6575=hyperGTest(params)

GO_6575 <- summary(overRepresented_6575)
write.csv(GO_6575, file="GO_6575")

# clustering/create heat plot
expmatrix_6575 <- exprs_6575[gene_list_6575_narm,]
expmatrix_cor_6575<-cor(t(expmatrix_6575))
expmatrix_dist_6575<-as.dist(1-expmatrix_cor_6575)

# cluster
expmatrix_hclust_6575<-hclust(
    expmatrix_dist_6575, method="ave"
)

expmatrix_hclust_groups_6575 = cutree(
    expmatrix_hclust_6575, k=2
)

names(
    which(expmatrix_hclust_groups_6575==1)
) -> group1_hclust_6575

names(
    which(expmatrix_hclust_groups_6575==2)
) -> group2_hclust_6575
pdf(file="boxplot1_6575.pdf")
boxplot(pair_6575@means[group1_hclust_6575,], main="Cluster 1: GSE6575",ylab="Mean expression",xlab="Group")
dev.off()
pdf(file="boxplot2_6575.pdf")
boxplot(pair_6575@means[group2_hclust_6575,], main="Cluster 2: GSE6575",ylab="Mean expression",xlab="Group")
dev.off()

sil_6575 <- silhouette(expmatrix_hclust_groups_6575, dist=expmatrix_dist_6575)
pdf("silhoutte_6575.pdf")
plot(sil_6575, col="black", main="Silhouette plot of GSE6575")
dev.off()


# create a heatmap
library(gplots)
hclust2 <- function(x, method="average", ...) {
    hclust(x, method=method, ...)
}
dist2 <- function(x, ...) {
    as.dist(1-cor(t(x), method="pearson"))
}
pdf("heatmap_6575.pdf")
heatmap.2(expmatrix_6575,
          col=redgreen(75),
          hclustfun=hclust2,
          distfun=dist2,
          scale="row",
          cexCol=0.6,
          Colv=TRUE,
          sepcolor="black",
          dendrogram="both",
          labCol=sapply(colnames(expmatrix_6575), function(x){covdesc_final_6575[covdesc_final_6575[,1]==x,2]}),
          key=TRUE,
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow = 0.4,
          main="Heatmap of GSE6575"
)
dev.off()

# get list of significant genes
geneSymbols_6575 <- unlist(mget(gene_list_6575_narm, hgu133plus2SYMBOL, ifnotfound = NA))
geneSymbols_6575_narm <- geneSymbols_6575[!is.na(geneSymbols_6575)]
merge(x=geneSymbols_6575_narm, y=pair_filt_6575@tt, by=union(names(x), names(y)))
pvals_sig <- vector()
i=0
for(aid in names(geneSymbols_6575_narm)){
    i = i+1
    indx <- match(aid, names(pair_filt_6575@tt))
    pvals_sig[i] <- pair_filt_6575@tt[indx]
}
pvals_sig <- as.numeric(pvals_sig)

geneSymbols_pval_6575 <- cbind(geneSymbols_6575_narm, as.numeric(pvals_sig))
geneSymbols_pval_6575 <- geneSymbols_pval_6575[order(geneSymbols_pval_6575[,2]),]
write.csv(file="significantGenes_6575.csv", geneSymbols_pval_6575)

# perform GO enrichment on significant genes
gene_list_6575_eid <- unlist(mget(gene_list_6575_narm, hgu133plus2ENTREZID, ifnotfound = NA))
gene_list_6575_eid_narm <- gene_list_6575_eid[!is.na(gene_list_6575_eid)]


params=new("GOHyperGParams", 
           geneIds=gene_list_6575_eid_narm,
           universeGeneIds=universe_eid_6575, 
           annotation="hgu133plus2", 
           ontology="BP",
           pvalueCutoff=0.001, 
           conditional=FALSE, 
           testDirection="over")

(overRepresented_6575=hyperGTest(params))
write.csv(summary(overRepresented_6575), file="goEnrichment_6575.csv")
# DONE WITH 6575




#########################################
# GSE7329 #

#getting GSE7329, "Gene expression profiles of lymphoblastoid cells (autism)"

# it appears that when you download GSE7316_RAW.tar from GEO, you actually get a bunch of txt.gz files, which are already normalized. 
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
geo_7329 <- getGEO('GSE7329', destdir=getwd()) # data is already log normalized based on cy3-cy5 ratios.
exprs_geo_7329 <- exprs(geo_7329[[1]])
file_designations_7329 = colnames(exprs_geo_7329)
design_7329 = c(rep('autism_dup15',7),rep('autism_fmr1',8),rep('control',15))
design_A_7329 = character()
for (i in 1:length(file_designations_7329)) {
    gsm_7329 = getGEO(file_designations_7329[i], destdir=getwd())
    design_A_7329[i] = strsplit(Meta(gsm_7329)$title, split=" ")[[1]][1]
}

# covdesc_7329 <- cbind(files_7329, design_7329)
# colnames(covdesc_7329) <- c("", "disease")
# write.table(x=covdesc_7329, file="covdesc_7329", sep="   ", row.names=F, quote=F)
# 
# Can't use read.affy, not a CEL format file.  File is already normalized.
# read.affy("covdesc_7329", compress=T, verbose=T)

# get log fold change
range_7329 <- apply(exprs_geo_7329, 1, getRange, design_7329)

# perform anova
anova_7329 <- apply(exprs_geo_7329, 1, doAnova, design_7329) #returns ~44,000 elements
fdr_7329 = p.adjust(anova_7329, method="BH") # correct for multiple hypothesis testing (FDR)

# check adjusted p-value - SJ
sig_genes_fdr_7329<- fdr_7329[fdr_7329 < 0.05]

length(sig_genes_fdr_7329) # approximately 3399

# get agilent ids for all id'd sig_genes
gene_names_7329_fdr <- names(sig_genes_fdr_7329)
gene_names_7329_narm <- gene_names_7329_fdr[!is.na(gene_names_7329_fdr)]

# get platform info
gpl1708 <- getGEO(filename='gpl1708.soft', destdir=getwd())

# get gene names - SJ
universe_7329 <- rownames(exprs_geo_7329)

# get entrezIDs from gene names
    # first get spot IDs from agilent
conversion_table_1708 <- Table(gpl1708)[,c(1,5,13,10)] # create conversion table
universe_SPOTID_7329<-conversion_table_1708[conversion_table_1708[,1] %in% universe_7329,2] # universe spot ids
universe_SPOTID_7329_narm <- as.character(universe_SPOTID_7329[!is.na(universe_SPOTID_7329)]) # remove NA
sig_matrix_7329 <- sig_genes_fdr_7329
names(sig_genes_fdr_7329) <- conversion_table_1708[conversion_table_1708[,1] %in% names(sig_genes_fdr_7329),2] # get spot ids for sig genes

SPOTID_7329_fdr <- names(sig_genes_fdr_7329[!is.na(names(sig_genes_fdr_7329))]) # remove NAs from spot ids

# get entrez IDs
universe_eid_7329 <- unlist(mget(universe_SPOTID_7329_narm, hgug4112aENTREZID, ifnotfound = NA)) # get entrez ids from db for universe
sig_genes_fdr_eid_7329 <- unlist(mget(SPOTID_7329_fdr, hgug4112aENTREZID, ifnotfound=NA)) # get entrez ids from db for sig genes

# go term analysis
params=new("GOHyperGParams", geneIds=sig_genes_fdr_eid_7329,
           universeGeneIds=universe_eid_7329, annotation="hgug4112a", ontology="BP",
           pvalueCutoff=0.001, conditional=FALSE, testDirection="over")

(overRepresented_7329=hyperGTest(params))

GO_7329 <- summary(overRepresented_7329)
write.csv(GO_7329, file="GO_7329.csv")

# clustering/create heat plot
expmatrix_7329 <- exprs_geo_7329[gene_names_7329_narm,]
colnames(expmatrix_7329) <- design_7329
expmatrix_cor_7329 <- cor(t(expmatrix_7329))
expmatrix_dist_7329<-as.dist(1-expmatrix_cor_7329)

# cluster
expmatrix_hclust_7329<-hclust(
    expmatrix_dist_7329, method="ave"
)

expmatrix_hclust_groups_7329 = cutree(
    expmatrix_hclust_7329, k=2
)

names(
    which(expmatrix_hclust_groups_7329==1)
) ->group1_hclust_7329

names(
    which(expmatrix_hclust_groups_7329==2)
) ->group2_hclust_7329


pdf("boxplot1_7329.pdf")
boxplot(exprs_geo_7329[group1_hclust_7329,], main="Cluster 1: GSE7329", ylab="Mean expression", xlab="Sample ID")
dev.off()
pdf("boxplot2_7329.pdf")
boxplot(exprs_geo_7329[group2_hclust_7329,], main="Cluster 2: GSE7329", ylab="Mean expression", xlab="Sample ID")
dev.off()

#silhoutte plot
sil_7329 <- silhouette(expmatrix_hclust_groups_7329, dist=expmatrix_dist_7329)
pdf("silhoutte_7329.pdf")
plot(sil_7329, col="black", main="Silhouette plot of GSE7329")
dev.off()

# get significant genes' symbol names based on cluster
clust1_symbols <- conversion_table_1708[conversion_table_1708[,1] %in% group1_hclust_7329, 4]
clust2_symbols <- conversion_table_1708[conversion_table_1708[,1] %in% group2_hclust_7329, 4]
names(clust1_symbols) <- group1_hclust_7329
names(clust2_symbols) <- group2_hclust_7329
clust1_symbols <- clust1_symbols[clust1_symbols!=""]
clust2_symbols <- clust2_symbols[clust2_symbols!=""]
length(clust1_symbols)
length(clust2_symbols)

# get significant genes' p-values
clust1_idx <- match(names(clust1_symbols), table=names(sig_matrix_7329))
clust2_idx <- match(names(clust2_symbols), table=names(sig_matrix_7329))
pvals_clust1 <- sig_matrix_7329[clust1_idx]
pvals_clust2 <- sig_matrix_7329[clust2_idx]

clust1_data <- cbind(as.character(clust1_symbols), pvals_clust1)
clust2_data <- cbind(as.character(clust2_symbols), pvals_clust2)

clust1_data <- clust1_data[order(clust1_data[,2]),]
clust2_data <- clust2_data[order(clust2_data[,2]),]

write.csv(clust1_data, file="repressedGenes_clust1_7329.csv")
write.csv(clust2_data, file="repressedGenes_clust2_7329.csv")

# create a heatmap
library(gplots)
hclust2 <- function(x, method="average", ...) {
    hclust(x, method=method, ...)
}
dist2 <- function(x, ...) {
    as.dist(1-cor(t(x), method="spearman"))
}
pdf("heatmap_7329.pdf")
heatmap.2(expmatrix_7329,
          col=redgreen(75),
          hclustfun=hclust2,
          distfun=dist2,
          scale="row",
          cexCol=0.65,
          Colv=TRUE,
          sepcolor="black",
          dendrogram="both",
          key=TRUE,
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow = 0.5,
          main="Heatmap of GSE7329"
)
dev.off()



##################
# get gse25507
# Download file, if not already downloaded
setwd(wd)
if(!file.exists("gse25507.tar") | !dir.exists("gse25507")){
    download.file("http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE25507&format=file", destfile="gse25507.tar")
    untar("gse25507.tar", exdir="gse25507")
}
setwd("gse25507")

# Get annotations from getGEO, get descriptions for each sample
geo25507 <- getGEO("GSE25507", destdir=getwd())
exprs_geo25507<- exprs(geo25507[[1]])
file_designations_25507 = colnames(exprs_geo25507)
design_25507 = vector(length=length(file_designations_25507))
for (i in 1:length(file_designations_25507)){
    filename = strsplit(file_designations_25507[i], split="[.]")[[1]][1]
    print(filename)
    gsm = getGEO(GEO=filename)
    design_25507[i] = strsplit(Meta(gsm)$characteristics_ch1[2], split=" ")[[1]][2]
}

# get file names
files_25507 <- sapply(file_designations_25507, function(x){ paste0(x,".CEL.gz")})

# create covdesc file
covdesc_25507 <- cbind(files_25507, design_25507)
colnames(covdesc_25507) <- c("", "disease")
write.table(x=covdesc_25507, file="covdesc_25507", sep="   ", row.names=F, quote=F)

raw_25507 <- read.affy("covdesc_25507", verbose = TRUE)
norm_gse25507 <- mas5(raw_25507)
exprs_gse25507 <- exprs(norm_gse25507)

# Perform pairwise comparison
pair_25507 <- pairwise.comparison(norm_gse25507, group="disease", members=c("autism", "control"), raw_25507)

pair_filt_25507<-pairwise.filter(pair_25507, min.present.no=3, present.by.group=T, fc=log2(1.5), tt=0.05)

gene_list_25507<-rownames(pair_filt_25507@means)
length(gene_list_25507) #5328
# go term enrichment
universe_25507 <- rownames(exprs_gse25507)
universe_25507_eid <- unlist(mget(universe_25507, hgu133plus2ENTREZID))
gene_eid_25507 <- unlist(mget(gene_list_25507[!is.na(gene_list_25507)], hgu133plus2ENTREZID))
gene_symbols_25507 <- unlist(mget(gene_list_25507[!is.na(gene_list_25507)], hgu133plus2SYMBOL))

params_25507=new("GOHyperGParams", geneIds=gene_eid_25507,
           universeGeneIds=universe_25507_eid, annotation="hgu133plus2", ontology="BP",
           pvalueCutoff=0.001, conditional=FALSE, testDirection="over")

(overRepresented_25507=hyperGTest(params_25597))

GO_25507 <- summary(overRepresented_25507)
write.csv(GO_25507, file="GO_25507.csv")
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
pdf("boxplot1_25507.pdf")
boxplot(pair_25507@means[group1_hclust_25507,], xlab="Group", ylab="Mean Expression", main="Cluster 1: GSE25507")
dev.off()
pdf("boxplot2_25507.pdf")
boxplot(pair_25507@means[group2_hclust_25507,], xlab="Group", ylab="Mean Expression", main="Cluster 2: GSE25507")
dev.off()
#silhoutte plot
sil_25507 <- silhouette(expmatrix_hclust_groups_25507, dist=expmatrix_dist_25507)
pdf("silhoutte_25507.pdf")
plot(sil_25507, col="black", main="Silhouette Plot of GSE25507")
dev.off()


# create a heatmap
library(gplots)
hclust2 <- function(x, method="average", ...) {
    hclust(x, method=method, ...)
}
dist2 <- function(x, ...) {
    as.dist(1-cor(t(x), method="pearson"))
}
pdf("heatmap_25507.pdf")
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
          cexRow = 0.4,
          main="Heatmap of GSE25507"
          
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

