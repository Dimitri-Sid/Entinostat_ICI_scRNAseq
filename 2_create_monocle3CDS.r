library(monocle3)
library(Matrix)
library(dplyr) 
#BiocManager::install("GEOquery")
library(GEOquery)

sessionInfo()

ConvertToMtx <- function(AnnData_File, Mtx_File) {
    suppressMessages(library(Matrix))
    suppressMessages(library(rhdf5))
    h5disableFileLocking()
    raw_data <- as.numeric(h5read(AnnData_File, "/X/data")) 
    raw_data <- as.numeric(format(raw_data, scientific = FALSE))
    indices <- as.numeric(h5read(AnnData_File, "/X/indices"))
    indptr <- as.integer(h5read(AnnData_File, "/X/indptr"))
    counts <- sparseMatrix(i=indices+1, p=indptr, x=raw_data)
    writeMM(counts, file=Mtx_File)
}

ConvertToMtx("aggr&pilot_filtered.h5ad", "aggr&pilot_filtered.mtx")

mtx <- readMM("aggr&pilot_filtered.mtx")

dim(mtx)

cell_metadata <- read.table("filtered_barcodes.txt")
gene_ids<-read.table("filtered_genes_IDs.txt")
gene_short_names<-read.table("filtered_genes_short_name.txt")

gene_annotation <- data.frame("full ID" = gene_ids,"gene_short_name" = gene_short_names)
colnames(gene_annotation) <- c("ensemblid", "gene_short_name")
rownames(gene_annotation)<-gene_annotation$ensemblid
rownames(cell_metadata)<-cell_metadata$V1

#create monocle3 cds object
cds <- new_cell_data_set(mtx,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
rm(mtx)
cds

cds <- estimate_size_factors(cds)

cds <- preprocess_cds(cds, num_dim = 25)

#unzip barcodes
gzipF <- list.files(path = "./data/aggr/aggr_sample_barcodes/", pattern = "*.gz", full.names = TRUE)
plyr::ldply(.data = gzipF, .fun = gunzip)

# Adding metadata by sample, 16 tumors from aggr
# modify to match pData in order to annotate barcodes by sample
AE3DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AE3DGELIB.tsv", header=FALSE)
AE3DGELIB$V1 <- gsub("-1", "", AE3DGELIB$V1)
AE3DGELIB$V1 <- paste(AE3DGELIB$V1,'-1-0', sep = "")
rownames(AE3DGELIB)<-AE3DGELIB$V1

AEC13DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AEC13DGELIB.tsv", header=FALSE)
AEC13DGELIB$V1 <- gsub("-1", "", AEC13DGELIB$V1)
AEC13DGELIB$V1 <- paste(AEC13DGELIB$V1,'-2-0', sep = "")
rownames(AEC13DGELIB)<-AEC13DGELIB$V1

AEC23DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AEC23DGELIB.tsv", header=FALSE)
AEC23DGELIB$V1 <- gsub("-1", "", AEC23DGELIB$V1)
AEC23DGELIB$V1 <- paste(AEC23DGELIB$V1,'-3-0', sep = "")
rownames(AEC23DGELIB)<-AEC23DGELIB$V1

AEP13DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AEP13DGELIB.tsv", header=FALSE)
AEP13DGELIB $V1 <- gsub("-1", "", AEP13DGELIB$V1)
AEP13DGELIB $V1 <- paste(AEP13DGELIB$V1,'-4-0', sep = "")
rownames(AEP13DGELIB)<-AEP13DGELIB$V1

AEP23DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AEP23DGELIB.tsv", header=FALSE)
AEP23DGELIB$V1 <- gsub("-1", "", AEP23DGELIB$V1)
AEP23DGELIB$V1 <- paste(AEP23DGELIB$V1,'-5-0', sep = "")
rownames(AEP23DGELIB)<-AEP23DGELIB$V1

AEPC13DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AEPC13DGELIB.tsv", header=FALSE)
AEPC13DGELIB$V1 <- gsub("-1", "", AEPC13DGELIB$V1)
AEPC13DGELIB$V1 <- paste(AEPC13DGELIB$V1,'-6-0', sep = "")
rownames(AEPC13DGELIB)<-AEPC13DGELIB$V1

AEPC23DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AEPC23DGELIB.tsv", header=FALSE)
AEPC23DGELIB $V1 <- gsub("-1", "", AEPC23DGELIB$V1)
AEPC23DGELIB $V1 <- paste(AEPC23DGELIB $V1,'-7-0', sep = "")
rownames(AEPC23DGELIB)<-AEPC23DGELIB$V1

AV3DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/AV3DGELIB.tsv", header=FALSE)
AV3DGELIB$V1 <- gsub("-1", "", AV3DGELIB$V1)
AV3DGELIB$V1 <- paste(AV3DGELIB$V1,'-8-0', sep = "")
rownames(AV3DGELIB)<-AV3DGELIB$V1

BE3DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BE3DGELIB.tsv", header=FALSE)
BE3DGELIB$V1 <- gsub("-1", "", BE3DGELIB$V1)
BE3DGELIB$V1 <- paste(BE3DGELIB$V1,'-9-0', sep = "")
rownames(BE3DGELIB)<-BE3DGELIB$V1

BEC13DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BEC13DGELIB.tsv", header=FALSE)
BEC13DGELIB$V1 <- gsub("-1", "", BEC13DGELIB$V1)
BEC13DGELIB$V1 <- paste(BEC13DGELIB$V1,'-10-0', sep = "")
rownames(BEC13DGELIB)<-BEC13DGELIB$V1

BEC23DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BEC23DGELIB.tsv", header=FALSE)
BEC23DGELIB$V1 <- gsub("-1", "", BEC23DGELIB$V1)
BEC23DGELIB$V1 <- paste(BEC23DGELIB$V1,'-11-0', sep = "")
rownames(BEC23DGELIB)<-BEC23DGELIB$V1

BEP13DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BEP13DGELIB.tsv", header=FALSE)
BEP13DGELIB$V1 <- gsub("-1", "", BEP13DGELIB$V1)
BEP13DGELIB$V1 <- paste(BEP13DGELIB$V1,'-12-0', sep = "")
rownames(BEP13DGELIB)<-BEP13DGELIB$V1

BEP23DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BEP23DGELIB.tsv", header=FALSE)
BEP23DGELIB$V1 <- gsub("-1", "", BEP23DGELIB$V1)
BEP23DGELIB$V1 <- paste(BEP23DGELIB$V1,'-13-0', sep = "")
rownames(BEP23DGELIB)<-BEP23DGELIB$V1

BEPC13DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BEPC13DGELIB.tsv", header=FALSE)
BEPC13DGELIB$V1 <- gsub("-1", "", BEPC13DGELIB$V1)
BEPC13DGELIB$V1 <- paste(BEPC13DGELIB$V1,'-14-0', sep = "")
rownames(BEPC13DGELIB)<-BEPC13DGELIB$V1

BEPC23DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BEPC23DGELIB.tsv", header=FALSE)
BEPC23DGELIB$V1 <- gsub("-1", "", BEPC23DGELIB$V1)
BEPC23DGELIB$V1 <- paste(BEPC23DGELIB$V1,'-15-0', sep = "")
rownames(BEPC23DGELIB)<-BEPC23DGELIB$V1

BV3DGELIB <- read.csv("./data/aggr/aggr_sample_barcodes/BV3DGELIB.tsv", header=FALSE)
BV3DGELIB$V1 <- gsub("-1", "", BV3DGELIB$V1)
BV3DGELIB$V1 <- paste(BV3DGELIB$V1,'-16-0', sep = "")
rownames(BV3DGELIB)<-BV3DGELIB$V1

#read in pilot barcodes, 4 more tumors (20 total, 5 per treatment group)
barcodes_wtile1<-read.table(file='./data/wtile1/barcodes.tsv')
barcodes_wtile1$V1 <- paste(barcodes_wtile1$V1,"-1", sep = "")
rownames(barcodes_wtile1)<-barcodes_wtile1$V1

barcodes_wtile2<-read.table(file='./data/wtile2/barcodes.tsv')
barcodes_wtile2$V1 <- paste(barcodes_wtile2$V1, "-2", sep = "")
rownames(barcodes_wtile2)<-barcodes_wtile2$V1

barcodes_wtilv1<-read.table(file='./data/wtilv1/barcodes.tsv')
barcodes_wtilv1$V1 <- paste(barcodes_wtilv1$V1, "-3", sep = "")
rownames(barcodes_wtilv1)<-barcodes_wtilv1$V1

barcodes_wtilv2<-read.table(file='./data/wtilv2/barcodes.tsv')
barcodes_wtilv2$V1 <- paste(barcodes_wtilv2$V1, "-4", sep = "")
rownames(barcodes_wtilv2)<-barcodes_wtilv2$V1

#add sample metadata
pData(cds)$sample <- ""
barcodes <- pData(cds)
barcodes$sample <- ""

barcodes[rownames(barcodes) %in% rownames(AV3DGELIB),]$sample <- "AV3DGELIB"
barcodes[rownames(barcodes) %in% rownames(BE3DGELIB),]$sample <- "BE3DGELIB"
barcodes[rownames(barcodes) %in% rownames(BEC13DGELIB),]$sample <- "BEC13DGELIB"
barcodes[rownames(barcodes) %in% rownames(BEC23DGELIB),]$sample <- "BEC23DGELIB"
barcodes[rownames(barcodes) %in% rownames(BEP13DGELIB),]$sample <- "BEP13DGELIB"
barcodes[rownames(barcodes) %in% rownames(BEP23DGELIB),]$sample <- "BEP23DGELIB"
barcodes[rownames(barcodes) %in% rownames(BEPC13DGELIB),]$sample <- "BEPC13DGELIB"
barcodes[rownames(barcodes) %in% rownames(BEPC23DGELIB),]$sample <- "BEPC23DGELIB"
barcodes[rownames(barcodes) %in% rownames(BV3DGELIB),]$sample <- "BV3DGELIB"
barcodes[rownames(barcodes) %in% rownames(AEPC13DGELIB),]$sample <- "AEPC13DGELIB"
barcodes[rownames(barcodes) %in% rownames(AE3DGELIB),]$sample <- "AE3DGELIB"
barcodes[rownames(barcodes) %in% rownames(AEC13DGELIB),]$sample <- "AEC13DGELIB"
barcodes[rownames(barcodes) %in% rownames(AEC23DGELIB),]$sample <- "AEC23DGELIB"
barcodes[rownames(barcodes) %in% rownames(AEP13DGELIB),]$sample <- "AEP13DGELIB"
barcodes[rownames(barcodes) %in% rownames(AEP23DGELIB),]$sample <- "AEP23DGELIB"
barcodes[rownames(barcodes) %in% rownames(AEPC23DGELIB),]$sample <- "AEPC23DGELIB"
barcodes[rownames(barcodes) %in% rownames(barcodes_wtile1),]$sample <- "barcodes_wtile1"
barcodes[rownames(barcodes) %in% rownames(barcodes_wtile2),]$sample <- "barcodes_wtile2"
barcodes[rownames(barcodes) %in% rownames(barcodes_wtilv1),]$sample <- "barcodes_wtilv1"
barcodes[rownames(barcodes) %in% rownames(barcodes_wtilv2),]$sample <- "barcodes_wtilv2"

pData(cds)$sample <- barcodes$sample

#add run annotation
pData(cds)$run <- ''
pData(cds)[pData(cds)$sample=="AE3DGELIB"|pData(cds)$sample=="AEC13DGELIB"|pData(cds)$sample=="AEC23DGELIB"|pData(cds)$sample=="AEP13DGELIB"|pData(cds)$sample=="AEP23DGELIB"|pData(cds)$sample=="AEPC13DGELIB"|pData(cds)$sample=="AEPC23DGELIB"|pData(cds)$sample=="AV3DGELIB",]$run <- "A"
pData(cds)[pData(cds)$sample=="BE3DGELIB" |pData(cds)$sample=="BEC13DGELIB"|pData(cds)$sample=="BEC23DGELIB"|pData(cds)$sample=="BEP13DGELIB"|pData(cds)$sample=="BEP23DGELIB"|pData(cds)$sample=="BEPC13DGELIB"|pData(cds)$sample=="BEPC23DGELIB"|pData(cds)$sample=="BV3DGELIB",]$run <- "B"

pData(cds)[pData(cds)$sample=="barcodes_wtile1",]$run <- "pilot1"
pData(cds)[pData(cds)$sample=="barcodes_wtile2",]$run <- "pilot2"
pData(cds)[pData(cds)$sample=="barcodes_wtilv1",]$run <- "pilot1"
pData(cds)[pData(cds)$sample=="barcodes_wtilv2",]$run <- "pilot2"

#add treatment annotation
pData(cds)$treatment <- ''
pData(cds)[pData(cds)$sample=="BV3DGELIB"|pData(cds)$sample=="AV3DGELIB"|pData(cds)$sample=="barcodes_wtilv1"|pData(cds)$sample=="barcodes_wtilv2",]$treatment <- "V"
pData(cds)[pData(cds)$sample=="AE3DGELIB"|pData(cds)$sample=="BE3DGELIB"|pData(cds)$sample=="barcodes_wtile1"|pData(cds)$sample=="barcodes_wtile2",]$treatment <- "E"
pData(cds)[pData(cds)$sample=="AEC23DGELIB"|pData(cds)$sample=="BEC13DGELIB"|pData(cds)$sample=="BEC23DGELIB"|pData(cds)$sample=="AEC13DGELIB",]$treatment <- "EC"
pData(cds)[pData(cds)$sample=="BEP13DGELIB"|pData(cds)$sample=="BEP23DGELIB"|pData(cds)$sample=="AEP13DGELIB"|pData(cds)$sample=="AEP23DGELIB",]$treatment <- "EP"
pData(cds)[pData(cds)$sample=="BEPC13DGELIB"|pData(cds)$sample=="BEPC23DGELIB"|pData(cds)$sample=="AEPC13DGELIB"|pData(cds)$sample=="AEPC23DGELIB",]$treatment <- "EPC"

#batch correct for "run" or library 
cds <- align_cds(cds, alignment_group = "run")
cds <- reduce_dimension(cds, max_components = 2, reduction_method = "UMAP") 
cds <- cluster_cells(cds, reduction_method = "UMAP", verbose = TRUE) 
plot_cells(cds, color_cells_by="cluster", group_cells_by="cluster")

# store clusters as metadata
pData(cds)$clusters <- ''
pData(cds)$clusters<-clusters(cds)
pData(cds)$partitions <- ''
pData(cds)$partitions<-partitions(cds)

plot_cells(cds, group_cells_by="partition", color_cells_by="treatment", cell_size = 0.3, label_cell_groups=FALSE,labels_per_group = 0, group_label_size = 10)

save(cds, file = "CDS.rda")


