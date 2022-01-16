# Supplemental and Pre-processing code

rm(list = ls())  # Clear the environment
options(warn=-1) # Turn off warning message globally
library(monocle3) # Load Monoc
library(dplyr)
library(igraph)
library(reticulate)
library(devtools)
library(reshape2)
library(ggplot2)
library(garnett)
library(pheatmap)
library(viridis)
library(gdata)
library(RColorBrewer)
load("CDS.rda") 

###########################################
# Broad cell types
###########################################
tcell_genes <- c("Cd8a","Cd8","Cd2","Cd3d","Cd3e","Cd3g","Cd28","Cd4","Cd8b1","Foxp3","Itgae")
b_cells <- c("Cd19","Cd79b","Cd45r","Cd79a","B220")
dc_genes <- c("Itgax", "Flt3", "Csf1r")
treg<-c("Cd3e", "Cd4","Foxp3","Itgae","Cd25","Il2r","Il-2r", "Il2ra", "Il-2ra") #must have CD4, FOXP3, and CD25 (il2ra)
nk_genes<- c("Nkg7","Klrb1", "Klrk1", "Ncr1", "Klrc1", "Cd56", "Fcgr3a", "Fcg3","Fcer1g")
t_active<-c("Hla","Cd38","Cd39","Cd69","Tcra","Tcrb","Tcr","H2-eb1","H2eb1")
t_exhausted <- c("Lag3","Cd244","Eomes","Ptger4")
macro <- c("F4/80","Adgre1","Cd68","Cd45","Ptprc","Cd64","Fcgr1","Mertk")
monocyte <- c("Cd14","Ccr2","Cd115","Csf1r","Gr-1","Cd11b","Vegf","Cx3cr1","Cd44","Cd31","Cd43","Ccl2")
M1 <- c("Ccl2","Cxcl1","Cd163","Wfdc17","Tnfa","Tnfsf1a","Tnf")
M2 <- c("Arg1","Ccl22","Ccl18","Stat6","Cd206", "Mrc1")
mdsc_genes <- c("Itgam", "Ly6g", "Ly6c1","Ly6c2")  #must be "Cd3e", "Adgre1" negative 
neutrophil_genes <-c("Itgam","S100a9","Cxcr4","Ly6g")
monocyte_macrophage_genes <- c("Adgre1", "Cd68", "Ptprc", "Fcgr1", "Mertk", "Itgam", "Csf1r", "Ccr2")
b_cells<- c("CD19","B220", "CD45R", "CD45r" )
nk_genes_2 <-c("Itga2", "Itgam", "Il2rb", "Klrb1", "Klrk1", "Ncr1", "Cd3e")
hsc_genes <- c("Emcn", "Kdr", "Kit", "Ly6a", "Cd48", "Slamf1")
cancer_genes_1 <-c("Lcn2", "Wfdc2", "S100a6", "Csn3", "Krt8", "Epcam", "Spp1","Cdh1")
cancer_genes_2 <- c("mt-Nd1", "mt-Nd5", "mt-Co1", "mt-Co3", "mt-Atp8", "mt-Cytb")
CAFs<-c("Itga3", "Itga5", "Oxtr", "Wnt5b", "Bcar1", "Fzd1")
VCAF <- c("Esam", "Gng11", "Higd1b", "Cox4i2", "Cygb", "Gja4", "Eng")
MCAF <- c("Dcn", "Col12a1", "Mmp2", "Lum", "Mrc2", "Bicc1", "Lrrc15", "Mfap5", "Col3A1", "Mmp14", "Spon1","Pdgfrl", "Serpinf1", "Lrp1", "Gfpt2", "Ctsk", "Cdh11", "Itgbl1", "Col6a2", "Postn", "Ccdc80", "Lox","Vcan", "Col1a1", "Fbn1", "Col1a2", "Pdpn", "Col6a1", "Fstl1", "Col5a2", "Aebp1")
DCAF <- c("Tspan2", "Reck")
mammary <- c("K6","K18","Krt19","Trp63","Pyn","Slp1","Ybx1","Eno1","CD49f", "Sca1", "Itgb1", "Krt14","Ly6a", "Cd14", "Krt8", "Krt6","Krt18")
g_mdsc <- c("Cd11b","Ly6g") #must be F40/80 negative
m_mdsc <- c("Cd11b") #must be Ly6G negative and F40/80 negative
t_memory <- c("Cd44","Cd62l")
mast <- c("Cd45","Cd117", "Fceri", "Itgb7")
plasma <- c("Cd138")
rbc <- c("Ter119, Ly76","Ly-76")
baso<- c("Cd200r3", "Fcer1a","Cd41","Cd49b") #must be Cd117 negative
eosino <- c("Cd11b","Cd193","F480","Singlec-f","Singlecf","Ccr3", "Il4ra") # must tbe MHCII negative
platelet_rest <- c("Cd41", "Cd9", "GpIa", "GpIIa", "Gp11b", "Gp11v", "Gpix", "Gpvi")
platelet_active <- c("Cd62p")
t_naive <- c("Cd45ra", "Ccr7", "Cd62l", "Cd127", "Cd132")
all_cell_type<-c("tcell_genes", "b_cells", "dc_genes", "treg", "nk_genes", "t_active", "t_exhausted", "macro","mono", 
                 "M1", "M2", "mdsc_genes", "neutrophil_genes", "monocyte_macrophage_genes", "b_cells", "nk_genes_2", 
                 "hsc_genes", "cancer_genes_1", "cancer_genes_2", "CAFs", "VCAF","MCAF", "DCAF", "mammary", "g_mdsc", "m_mdsc", 
                 "t_memory", "mast", "plasma", "rbc", "baso", "eosino", "platelet_rest", "platelet_active", "t_naive")
all_genes_grouped<-list(tcell_genes, b_cells, dc_genes, treg, nk_genes, t_active, t_exhausted, macro,monocyte, M1, M2, mdsc_genes, neutrophil_genes, monocyte_macrophage_genes, b_cells, nk_genes_2, hsc_genes, cancer_genes_1, cancer_genes_2, CAFs, VCAF, MCAF, DCAF, mammary, g_mdsc, m_mdsc, t_memory, mast, plasma, rbc, baso, eosino, platelet_rest, platelet_active, t_naive)

pdf("all_umaps.pdf")
for (i in seq(1,length(all_cell_type))) {   
  print(plot_cells(cds,genes = all_genes_grouped[[i]], label_cell_groups = F, cell_size = 0.5)+ggtitle(all_cell_type[i]))
} 
dev.off() 

# Assign cell types
colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                "1"="Cancer 1",
                                                "2"="Cancer 2",
                                                "3"="MDSCs",
                                                "4"="Macrophages/Monocytes",
                                                "5"="CAFs",
                                                "6"="Lymphoids",
                                                "7"="Doublets",
                                                "8"="Doublets",
)

# Remove doublets and cancer 2: 
cds <- cds[,pData(cds)$assigned_cell_type != "Cancer 2" | pData(cds)$assigned_cell_type != "doublets"]

save(cds,file="broad_cell_types_CDS.rda")

###########################################
# Immune populations - Myeloid and Lymphoid
###########################################
cds_myeloids <- cds[,pData(cds)$assigned_cell_type == "MDSCs" | pData(cds)$assigned_cell_type == "Monocytes/Macrophages"]
cds_lymphoids <- cds[,pData(cds)$assigned_cell_type == "Lymphoid"]
cds_myeloids <- reduce_dimension(cds_myeloids, max_components=2, reduction_method = c("UMAP"))
cds_lymphoids <- reduce_dimension(cds_lymphoids, max_components=2, reduction_method = c("UMAP"))
cds_myeloids = cluster_cells(cds_myeloids, resolution=1e-4)
cds_lymphoids = cluster_cells(cds_lymphoids, resolution=1e-3)

macrophage <- c("Il6","Mertk","Adgre1","Arg1","Il6","Tnf","Il1b","Stat3","Tgfb1","Tnf","Itgam","Gatm","C1qc","C1qb","C1qa","Ctss","Tmsb10","Cxcl2","Cd14","Cd2","Apoe","Lgmn","Lyz2","Cd34")
m1_vs_m2 <- c("Arg1","Il6","Tnf","Il1b")
tcell_genes <- c("Cd3e","Cd3d","Cd8a","Cd4","Foxp3","Itgae","Cd25","Il2ra","Ctla4")
b_cells <- c("Cd19","Cd79b","Cd45r","Cd79a","B220")
treg<-c("Cd3e","Tgfb1","CD127", "Cd4","Foxp3","Itgae","Il10","Cd25","Il2r","Il-2r", "Il2ra", "Il-2ra") #must have CD4, FOXP3, and CD25 (il2ra)
nkt_genes<- c("Nkg7","Klrb1", "Klrk1", "Ncr1", "Klrc1", "Cd56", "Fcgr3a", "Fcg3","Fcer1g")
t_active<-c("Hla","H2-Eb1","Erb1","Cd38","Cd39","Cd69","Tcra","Tcrb","Tcr","H2-eb1","H2eb1")
t_exhausted <- c("Lag3","Cd244","Eomes","Ptger4")
t_naive <- c("Cd45ra", "Ccr7", "Cd62l", "Cd127", "Cd132")
MDSC <- c("Ptges2","Stat3","Cd3d","Cd3e","Cd14","Cd15","Tnf","Cd33","Cd66b","S100a8","S100a9","Fth1","Cd34","Itgam","Ly6g","Ly6c2","Ly6c1","Adgre1") #appear different
MDSC_GvM <- c("Itga4","Ptges2","Cd3e","Cd14","Lyc2","S100a8","S100a9","Tnf","Itgam","Adgre1","Ly6g","Ly6c1","Ly6c2","Stat1","Stat2","Stat3","Arg1")
murine_mdsc<-c("Itgam","Adgre1","Ly6g","Ly6c1")
DC <- c("Flt3","Btla","Itgax")

pdf("macrophage_markers_in_myeloid.pdf")
plot_cells(cds_myeloids,genes = macrophage, label_cell_groups = F, cell_size = 0.5)+ggtitle("Markers")
dev.off()

pdf("t_cell_markers_in_lymphoid.pdf")
plot_cells(cds_lymphoids,genes = c(tcell_genes,t_reg,t_active,t_exhausted,t_naive), label_cell_groups = F, cell_size = 0.5)+ggtitle("Markers")
dev.off()

pdf("b_cell_markers_in_lymphoid.pdf")
plot_cells(cds_lymphoids,genes = b_cells, label_cell_groups = F, cell_size = 0.5)+ggtitle("Markers")
dev.off()

pdf("nkt_cell_markers_in_lymphoid.pdf")
plot_cells(cds_lymphoids,genes =nkt_genes, label_cell_groups = F, cell_size = 0.5)+ggtitle("Markers")
dev.off()

pdf("MDSC_markers_in_myeloid.pdf")
plot_cells(cds_myeloids,genes =MDSC, label_cell_groups = F, cell_size = 0.5)+ggtitle("Markers")
dev.off()

pdf("MDSC_markers_in_myeloid_mdsc.pdf")
plot_cells(cds_myeloids,genes =MDSC_GvM, label_cell_groups = F, cell_size = 0.8)+ggtitle("Markers")
dev.off()

#Rename clusters for myeloids
colData(cds_myeloids)$clusters <-  cds_myeloids@clusters$UMAP$clusters
colData(cds_myeloids)$clusters = dplyr::recode(colData(cds_myeloids)$clusters,
                                               "1"="Macrophages",
                                               "2"="G-MDSC",
                                               "3"="M-MDSC",
                                               "4"="Macrophages",
                                               "5"="Monocytes",
                                               "6"="Macrophages",
                                               "7"="Macrophages",
                                               "8"="DC",
)

#Rename clusters for lymphoids
colData(cds_lymphoids)$clusters <-  cds_lymphoids@clusters$UMAP$clusters
colData(cds_lymphoids)$clusters = dplyr::recode(colData(cds_lymphoids)$clusters,
                                                "1"="Treg",
                                                "2"="Th",
                                                "3"="Treg",
                                                "4"="B",
)
# Manually select the Tc cells 
cds_lymphoids_Tc <- choose_cells(cds_lymphoids) # to select Tc
colData(cds_lymphoids)[rownames(pData(cds_lymphoids_Tc)),"clusters"] <- as.factor("Tc")
levels(pData(cds_lymphoids)$clusters) <- c("Treg","Th","B","Tc")

# Select myeloid subpopulations
cds_myeloids_mdsc <- cds_myeloids[,pData(cds_myeloids)$assigned_cell_type == "MDSCs"]
cds_myeloids_macro <- cds_myeloids[,pData(cds_myeloids)$assigned_cell_type == "Monocytes/Macrophages"]
cds_myeloids_cDC <- cds_myeloids[,pData(cds_myeloids)$clusters == "DC"]
cds_myeloids_mdsc<-cluster_cells(cds_myeloids_mdsc)
cds_myeloids_macro<-cluster_cells(cds_myeloids_macro)
cds_myeloids_cDC<-cluster_cells(cds_myeloids_cDC)

# Calculate pseudotime
cds_lymphoids <- learn_graph(cds_lymphoids, use_partition = FALSE)
cds_myeloids <- learn_graph(cds_myeloids, use_partition = FALSE)
cds_myeloids_cDC <- learn_graph(cds_myeloids_cDC, use_partition = FALSE)
cds_myeloids_macro <- learn_graph(cds_myeloids_macro, use_partition = FALSE)
cds_myeloids_mdsc <- learn_graph(cds_myeloids_mdsc, use_partition = FALSE)

#save myeloid and lymphoid CDS objects
save(cds_myeloids,file="cds_myeloids.rda")
save(cds_lymphoids,file="cds_lymphoids.rda")