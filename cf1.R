setwd("/Users/menglinghe/Library/Mobile Documents/com~apple~CloudDocs/UTK/GRA-UTK/UTK_R/CF_project")

# load libraries
library(dplyr)
library(Seurat)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(harmony)
set.seed(412)

# Load the CF rds -----------------
cf <- readRDS("./results/cf.rds")


# data processing -----------------
cf.1 <- subset(cf, subset = nFeature_RNA>200 & nFeature_RNA < 7000 & percent.mt <5)
table(cf.1@meta.data$cell.type.ident)

# split object
cf.1[["RNA"]] <- split(cf.1[["RNA"]], f = cf.1$disease.ident)

# normalize data
cf.1 <- NormalizeData(cf.1)
# feature selection(find highly variable features)
cf.1 <- FindVariableFeatures(cf.1, selection.method = "vst", nfeatures = 3000)
# scale to mean 0 and variance 1
all.genes <- rownames(cf.1)
cf.1<-ScaleData(cf.1, features=all.genes, verbose=F)
# regress out percentage number of mitochondrial genes
cf.1 <-ScaleData(cf.1, vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)


# Perform Linear dimensionality reduction --------------
cf.1 <- RunPCA(cf.1, features = VariableFeatures(object = cf.1))
# visualize reduced dimensions
#VizDimLoadings(CF_integration.1, dims = 1:10, reduction = "pca")
# dimension heat map(first dimension- primary sources for heterogen.)
#cells and features ordered by PCA scores
#DimHeatmap(CF_integration.12, dims = 1:20, cells = 500, balanced = TRUE)
# Determine "dimensionality" of dataset (paper used 12)
ElbowPlot(cf.1, ndims=20, reduction="pca")

#  unintegrated clustering------------
cf.1 <- FindNeighbors(cf.1, dims = 1:15, reduction = "pca")
cf.1 <- FindClusters(cf.1, resolution = 0.5, cluster.name = "unintegrated_clusters")
cf.1 <- RunUMAP(cf.1, dims = 1:15, reduction = "pca", reduction.name = "umap.unintegrated")
cf.1 <- RunTSNE(cf.1, dims = 1:15, reduction="pca", reduction.name = "tsne.unintegrated")
#plot2 <- TSNEPlot(cf.1, reductions="tsne.unintegrated"),group.by = c("disease.ident", "unintegrated_clusters","disease.ident"))
plot1 <- DimPlot(cf.1, reduction = "umap.unintegrated", group.by = c("disease.ident", "cell.type.ident", "unintegrated_clusters"))

saveRDS(cf.1, file = "./results/cf_1.rds")



# Perform streamlined (one-line) integrative analysis--------------------
# try 3 methods
CF_integration.1 <- IntegrateLayers(
  object = cf.1, method = CCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.cca',
  verbose = FALSE)
CF_integration.1 <- IntegrateLayers(
  object = CF_integration.1, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.rpca',
  verbose = FALSE)
CF_integration.1 <- IntegrateLayers(
  object = CF_integration.1, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# CF_integration.1 <- IntegrateLayers(
#   object = CF_integration.1, method = FastMNNIntegration,
#   new.reduction = "integrated.mnn",
#   verbose = FALSE
# )
# CF_integration.1 <- IntegrateLayers(
#   object = CF_integration.1, method = scVIIntegration,
#   new.reduction = "integrated.scvi",
#   conda_env = "../miniconda3/envs/scvi-env", verbose = FALSE
# )


# clustering (3 methods)

CF_integration.1 <- FindNeighbors(CF_integration.1, reduction = "integrated.rpca", dims = 1:15)
CF_integration.1 <- FindClusters(CF_integration.1, resolution =c(0.3), cluster.name = c("rpca_clusters_0.3"))
CF_integration.1 <- RunUMAP(CF_integration.1, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
p2 <- DimPlot(
  CF_integration.1,
  reduction = "umap.rpca",
  group.by = c("disease.ident", "cell.type.ident"),
  combine = FALSE, label.size = 2,label = T
)


CF_integration.1 <- FindNeighbors(CF_integration.1, reduction = "integrated.cca", dims = 1:15)
CF_integration.1 <- FindClusters(CF_integration.1, resolution = 0.5,cluster.name = "cca_clusters")
CF_integration.1 <- RunUMAP(CF_integration.1,dims = 1:15,reduction = 'integrated.cca',reduction.name = "umap.cca")
p1 <- DimPlot(
  CF_integration.1,
  reduction = "umap.cca",
  group.by = c("disease.ident", "cell.type.ident"),
  combine = FALSE, label.size = 2,label = T
)

CF_integration.1 <- FindNeighbors(CF_integration.1, reduction = "integrated.rpca", dims = 1:15)
CF_integration.1 <- FindClusters(CF_integration.1, resolution =c(0.1,0.3,0.5), cluster.name = c("rpca_clusters_0.1","rpca_clusters_0.3","rpca_clusters_0.5"))
CF_integration.1 <- RunUMAP(CF_integration.1, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
p2 <- DimPlot(
  CF_integration.1,
  reduction = "umap.rpca",
  group.by = c("disease.ident", "cell.type.ident"),
  combine = FALSE, label.size = 2,label = T
)

CF_integration.1 <- FindNeighbors(CF_integration.1, reduction = "harmony", dims = 1:15)
CF_integration.1 <- FindClusters(CF_integration.1, resolution = 0.5, cluster.name = "harmony_clusters")
CF_integration.1 <- RunUMAP(CF_integration.1, reduction = "harmony", dims = 1:15, reduction.name = "umap.harmony")
p3 <- DimPlot(
  CF_integration.1,
  reduction = "harmony",
  group.by = c("disease.ident", "cell.type.ident"),
  combine = FALSE, label.size = 2,label = T
)

wrap_plots(c(p1, p2, p3), ncol = 2, byrow = T)
# choose RPCA


# p4 <- DimPlot(CF_integration.1, reduction = "umap.unintegrated", group.by = c("cca_clusters"))
# p5 <- DimPlot(CF_integration.1, reduction = "umap.cca", group.by = c("cca_clusters"))
# p6 <- DimPlot(CF_integration.1, reduction = "umap.rpca", group.by = c("cca_clusters"))
# p7 <- DimPlot(CF_integration.1, reduction = "umap.harmony", group.by = c("cca_clusters"))
# grid.arrange(p4,p5,p6,p7, ncol = 2, nrow = 2)

# show RPCA results in multiple reduction dim
p4_1 <- DimPlot(CF_integration.1, reduction = "umap.unintegrated", group.by = c("rpca_clusters_0.5"),label = T)
p5_1 <- DimPlot(CF_integration.1, reduction = "umap.cca", group.by = c("rpca_clusters_0.5"),label = T)
p6_1 <- DimPlot(CF_integration.1, reduction = "umap.rpca", group.by = c("rpca_clusters_0.5"),label = T)
p7_1 <- DimPlot(CF_integration.1, reduction = "umap.harmony", group.by = c("rpca_clusters_0.5"),label = T)
grid.arrange(p4_1,p5_1,p6_1,p7_1, ncol = 2, nrow = 2)


# compare integration result vs. cell.type 
plot_clustering_rpca <- DimPlot(CF_integration.1, reduction = "umap.rpca"
        , group.by = c("cell.type.ident","rpca_clusters_0.1"
                       ,"rpca_clusters_0.3","rpca_clusters_0.5")
        ,label = T)

plot_clustering_un <- DimPlot(CF_integration.1, reduction = "umap.unintegrated"
                                , group.by = c("cell.type.ident","rpca_clusters_0.1"
                                               ,"rpca_clusters_0.3","rpca_clusters_0.5")
                                ,label = T)


# # add resolution =0.2 to rpca
# CF_integration.1 <- FindNeighbors(CF_integration.1, reduction = "integrated.rpca", dims = 1:15)
# CF_integration.1 <- FindClusters(CF_integration.1, resolution =c(0.2), cluster.name = c("rpca_clusters_0.2"))
# CF_integration.1 <- RunUMAP(CF_integration.1, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")

# 
# 
# DimPlot(CF_integration.1, reduction = "umap.rpca"
#         , group.by = c("rpca_clusters_0.1","rpca_clusters_0.2"
#                        ,"rpca_clusters_0.3","rpca_clusters_0.5")
#         ,label = T)
# 















# table showing each cell by cluster ------------------------
CF_integration.1[["cell_bycluster"]] <- Idents(CF_integration.1)

CF_summary <- CF_integration.1@meta.data[c("disease.ident","rpca_clusters_0.1")]
tab <- table(CF_summary)
tab <- cbind(tab, Total = rowSums(tab))



# cell by cluster ------------------------
cluster_bycell <- CF_integration.1@meta.data["rpca_clusters_0.1"]

write.csv(cluster_bycell, "./results/cluster_bycell.csv", row.names=TRUE)



# findConserved markers (resolution=0.1)--------------------
CF_integration.1 <- JoinLayers(CF_integration.1)
saveRDS(CF_integration.1, file = "./results/CF_integration_1.rds")

CF_integration.1 <- readRDS("./results/CF_integration_1.rds")

Idents(CF_integration.1) <- CF_integration.1@meta.data$rpca_clusters_0.1

# CD68, SiglecF, CD163, CD206, EMR1, and CD11b
VlnPlot(CF_integration.1, features=c("CD68", "SiglecF","CD163"
                                     ,"CD206", "EMR1", "CD11b"))

# Macrophages markers
VlnPlot(CF_integration.1, features=c("FNIP2", "ALCAM"))


# plot 
DimPlot(CF_integration.1, reduction = "umap.unintegrated", split.by = "disease.ident")
DimPlot(CF_integration.1, reduction = "umap.rpca", split.by = "disease.ident")

# MRC1 or MARCO typically detected on alvMΦs 
VlnPlot(CF_integration.1, features=c("MRC1", "MARCO"))
# rename cluster 5 ident
CF_integration.1 <- RenameIdents(CF_integration.1, `5` = 'AM1')

# cluster0
marker_r1_0 <- FindConservedMarkers(CF_integration.1, ident.1 = 0, grouping.var = "disease.ident", verbose = FALSE)
head(marker_r1_0,20)
#FCGR3B,CXCR2,ALPL -- PMN by paper
VlnPlot(CF_integration.1, features=c("CXCR4", "IGF2R","ALPL"))
VlnPlot(CF_integration.1, features=c("FCGR3B", "CXCR2","ALPL"))
CF_integration.1 <- RenameIdents(CF_integration.1, `0` = 'PMN')


# cluster1&2 (both have Macrophages marker genes)
marker_r1_1 <- FindConservedMarkers(CF_integration.1, ident.1 = 1, grouping.var = "disease.ident", verbose = FALSE)
head(marker_r1_1)
# ANO5-- Macrophages;FAM89A--Macrophages;ADTRP--Macrophages
VlnPlot(CF_integration.1, features=c("MERTK", "CD64","SiglecF"))
FeaturePlot(CF_integration.1, features = c('MERTK')
            , reduction =  "umap.rpca", min.cutoff = 'q10')#cluster 2

# cluster2
marker_r1_2 <- FindConservedMarkers(CF_integration.1, ident.1 = 2, grouping.var = "disease.ident", verbose = FALSE)
head(marker_r1_2)
#GPR183-- Macrophages; TMSB4X--Macrophages

# most CF airway MΦs originate from recruited monocytes, whereas the majority of HC airway cells were bona fide tissue-resident alvMΦs
DimPlot(CF_integration.1,reduction = "umap.rpca"
        ,group.by = c("disease.ident", "rpca_clusters_0.1")
        , label.size = 4,label = T)

# These RLPs were defined by a high expression of mononuclear phagocyte–associated genes
FeaturePlot(CF_integration.1, features = c('LYZ','CTSB','CTSH','CTSL')
            ,reduction =  "umap.rpca", min.cutoff = 'q10')
FeaturePlot(CF_integration.1, features = c('LYZ','CTSB','CTSH','CTSL')
            ,reduction =  "umap.rpca", min.cutoff = 'q10')

CF_integration.1 <- RenameIdents(CF_integration.1, `1` = 'AM2')
CF_integration.1 <- RenameIdents(CF_integration.1, `2` = 'M/Momo')# 


# cluster3
marker_r1_3 <- FindConservedMarkers(CF_integration.1, ident.1 = 3, grouping.var = "disease.ident", verbose = FALSE)
head(marker_r1_3)
#KRT19
CF_integration.1 <- RenameIdents(CF_integration.1, `3` = 'epithelial')
VlnPlot(CF_integration.1, features=c("KRT19"))

# cluster4
marker_r1_4 <- FindConservedMarkers(CF_integration.1, ident.1 = 4, grouping.var = "disease.ident", verbose = FALSE)
head(marker_r1_4)
# PPP1R16B--Tcells; CARD11-- T/B;IKZF3--T/B; SLC38A1--T/B
CF_integration.1 <- RenameIdents(CF_integration.1, `4` = 'T/B')
VlnPlot(CF_integration.1, features=c("PPP1R16B","CARD11","SLC38A1"))
VlnPlot(CF_integration.1, features=c("PPP1R16B","CARD11","IKZF3"))

# cluster5
#marker_r1_5 <- FindConservedMarkers(CF_integration.1, ident.1 = 5, grouping.var = "disease.ident", verbose = FALSE)
marker_r1_5 <- FindMarkers(CF_integration.1, ident.1 = 'AM1', 
                           thresh.use = 0.25,  min.pct = 0.1,
                           only.pos = FALSE)
head(marker_r1_5)

# visualize the renamed cluster in unintegrated and rpca reduction dimension
plot_clustering1 <- DimPlot(CF_integration.1,reduction = c("umap.unintegrated")
        , label.size = 4,label = T)
plot_clustering2 <- DimPlot(CF_integration.1,reduction = c("umap.rpca")
          , label.size = 4,label = T)
plot_clustering1|plot_clustering2

plot_clustering1 <- DimPlot(CF_integration.1,reduction = c("umap.unintegrated")
                            ,label.size = 4,label = T,split.by = "disease.ident")
plot_clustering2 <- DimPlot(CF_integration.1,reduction = c("umap.rpca")
                            , label.size = 4,label = T,split.by = "disease.ident")

# get the first 20 markers of each cluster
N = 20
list_of_markers <- list(marker_r1_0, marker_r1_1, marker_r1_2
                        ,marker_r1_3,marker_r1_4,marker_r1_5)
list_of_markers <- lapply(list_of_markers, function(df) rownames(df)[1:N])
markers <- unlist(list_of_markers)
cluster <- rep(0:5, each = N)
# Combine into a DataFrame
marker.df <- data.frame(markers = markers, cluster = cluster)
write.csv(marker.df, "./results/genemarker_updated.csv", row.names=FALSE)






# findMarkers between conditions ---------------------
CF_integration.1$celltype.cnd <- paste0(CF_integration.1$rpca_clusters_0.1,'_', CF_integration.1$disease.ident)
View(CF_integration.1@meta.data)
Idents(CF_integration.1) <- CF_integration.1$celltype.cnd
DimPlot(CF_integration.1, reduction = 'umap.rpca', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(CF_integration.1, ident.1 = '1_CF', ident.2 = '1_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
head(marker_r1_1)
FeaturePlot(CF_integration.1, features = c('G0S2', 'PRKCB', 'ANO5','AC026369.3'), split.by = 'disease.ident', min.cutoff = 'q10')


















# save plots in PDF--------------------
#open PDF
pdf(file="./results/myplots.pdf")

#save plots to PDF
plot1
#grid.arrange(p4_1,p5_1,p6_1,p7_1, ncol = 2, nrow = 2)
plot_clustering_rpca
plot_clustering_un
plot_clustering1|plot_clustering2


#turn off PDF plotting
dev.off() 
