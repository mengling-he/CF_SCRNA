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

# load gene matrix
mat<-Seurat::ReadMtx("./CF_data/GSE145360_Sputum.data.processed.mtx",cells="./CF_data/GSE145360_Sputum.processed.barcodes.tsv",
                     features="./CF_data/GSE145360_Sputum.processed.genes.tsv", feature.column = 1)
#load metadata
meta <- data.frame(fread("https://cells.ucsc.edu/human-sputum/meta.tsv"), row.names=1)

# Initialize the Seurat object with the raw (non-normalized data).
# New seurat object for paper analysis
cf <- CreateSeuratObject(counts = mat, project = "cf", min.cells = 3, min.features = 200,
                           meta.data=meta)

Idents(cf) <- cf@meta.data$disease.ident

#table(cf@meta.data$cell.type.ident)

# 1. QC -------

# get percent mitochondrial genes
cf[["percent.mt"]] <- PercentageFeatureSet(cf, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(cf, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")
        #, ncol = 3
        ,group.by = NULL
        )

plot1 <- FeatureScatter(cf, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(cf, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")+geom_smooth(method = 'lm')
plot1 + plot2

saveRDS(cf, file = "./results/cf.rds")













# 2. Filtering -----------------
#subset data cf.1 and cf.2 is a subset of cf, but different filtering criteria
cf.1 <- subset(cf, subset = nFeature_RNA>200 & nFeature_RNA < 7000 & percent.mt <5)
saveRDS(cf.1, file = "../output/cf.1.rds")
cf.2 <- subset(cf, subset = nFeature_RNA>200 & nFeature_RNA < 2250 & percent.mt <20)
cell <- cf@meta.data$cell.type.ident
# compare the filtering ratio 
cell1 <- cf.1@meta.data$cell.type.ident
cell2 <- cf.2@meta.data$cell.type.ident
cell_compare12 <- data.frame(table(cell),table(cell1),table(cell2))
cell_compare12$ratio1 <- cell_compare12$Freq.1/cell_compare12$Freq
cell_compare12$ratio2 <- cell_compare12$Freq.2/cell_compare12$Freq

# 3. Normalize data ----------
cf.2 <- NormalizeData(cf.2)

# 4. Identify highly variable features --------------
cf.2 <- FindVariableFeatures(cf.2, selection.method = "vst", nfeatures = 3000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(cf.2), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(cf.2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

# 5. Scaling -------------
all.genes <- rownames(cf.2)
cf.2<-ScaleData(cf.2, features=all.genes, verbose=F)
#regress out percentage number of mitochondrial genes
cf.2 <-ScaleData(cf.2, vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)

#PCA dimensionaily reduction on subset of scaled variable features
cf.2 <- RunPCA(cf.2, features = VariableFeatures(object = cf.2))

#visualize reduced dimensions
VizDimLoadings(cf.2, dims = 1:10, reduction = "pca")

#dimension heat map(first dimension- primary sources for heterogen.)
#cells and features ordered by PCA scores
DimHeatmap(cf.2, dims = 1:20, cells = 500, balanced = TRUE)

#Determine "dimensionality" of dataset (paper used 12)
ElbowPlot(cf.2, ndims=20, reduction="pca")

#######Cluster the cells
#nearest neighbors
cf.2 <- FindNeighbors(cf.2, dims = 1:13)

#cluster with Louvain (no resolution size outlined in methods)
#cf.2 <- FindClusters(cf.2, resolution = 0.5)
cf.2 <- FindClusters(cf.2, resolution = c(0.5,0.7,1))

DimPlot(cf.2,group.by = "RNA_snn_res.0.5",label=TRUE)
DimPlot(cf.2,group.by = "RNA_snn_res.0.7",label=TRUE)
DimPlot(cf.2,group.by = "RNA_snn_res.1",label=TRUE)
# # setting identity of clusters based on dimplot of different resolution
Idents(cf.2) <- "RNA_snn_res.1"

##########Run non-linear dimensional reduction(UMAP/tSNE)
#try to place similar cells together in a low-dimensional space
#suggested use fo same PC's as clustering analysis

#UMAP-Unif.ManifoldAppx&projection
cf.2 <- RunUMAP(cf.2, dims = 1:12)

cf.2<-RunTSNE(cf.2, reduction="pca", dims = 1:12)

#visualize UMAP/tsne
DimPlot(cf.2, reduction = "umap", label=T)
TSNEPlot(cf.2, reduction="tsne", label=T)

before_dim <- DimPlot(cf.2, reduction = 'umap', group.by = 'disease.ident')
before_tsne <- DimPlot(cf.2, reduction = 'umap', group.by = 'disease.ident')


# 
# cf.3 <- RunUMAP(cf.2, dims = 1:12, reduction = "pca", reduction.name = "umap.unintegrated")
# DimPlot(cf.3, reduction = "umap.unintegrated", group.by = c("disease.ident", "seurat_clusters"))
# DimPlot(cf.3,reduction = "umap", label=T)#same result with cf.2

# #########Get cluster markers
# #find only positive markers with over 25% of cells in cluster expressing
# cf2.markers <- FindAllMarkers(cf.2, only.pos = TRUE, min.pct = 0.25,
#                               logfc.threshold = 0.25)
# CF2.MARKERS<-as.data.frame(cf2.markers %>%
#                              group_by(cluster) %>%
#                              slice_max(n = 20, order_by = avg_log2FC)
# )



# 
# 
# write.csv(CF2.MARKERS, file = "./results/CF2.MARKERS.csv", row.names = FALSE)
# 
# CF2.MARKERS_50<-as.data.frame(cf2.markers %>%
#                              group_by(cluster) %>%
#                              slice_max(n = 50, order_by = avg_log2FC)
# )
# write.csv(CF2.MARKERS_50, file = "./results/CF2.MARKERS.50.csv", row.names = FALSE)
# 










##############Integration
p1 <- DimPlot(cf.2,reduction='umap',group.by ='disease.ident' )
p2 <- DimPlot(cf.2,reduction='umap',group.by ='cell.type.ident' )
#grid.arrange(p1,p2,ncol = 2, nrow = 2)
p3 <- TSNEPlot(cf.2,reduction='tsne',group.by ='disease.ident' )
p4 <- TSNEPlot(cf.2,reduction='tsne',group.by ='cell.type.ident' )
grid.arrange(p1,p2,p3,p4, ncol = 2, nrow = 2)



cf.list <- SplitObject(cf.2, split.by = "disease.ident")

cf.list <- lapply(X = cf.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select integration features
features <- SelectIntegrationFeatures(object.list = cf.list)

immune.anchors <- FindIntegrationAnchors(object.list = cf.list, anchor.features = features)
immune.combined <- IntegrateData(anchorset = immune.anchors)

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <-ScaleData(immune.combined, vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined, ndims=20, reduction="pca")

# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:12)
immune.combined <- RunTSNE(immune.combined, reduction = "pca", dims = 1:12)

immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindClusters(immune.combined, resolution =  c(0.5,0.7,1))

# Visualization
p11 <- DimPlot(immune.combined, reduction = "umap", group.by = "disease.ident")
p21 <- DimPlot(immune.combined, reduction = "umap", group.by = "cell.type.ident")
p31 <- TSNEPlot(immune.combined, reduction = "tsne", group.by = "disease.ident")
p41 <- TSNEPlot(immune.combined, reduction = "tsne", group.by = "cell.type.ident")
grid.arrange(p11,p21,p31,p41, ncol = 2, nrow = 2)


p5 <- DimPlot(immune.combined, reduction = "umap", split.by = "disease.ident",label = TRUE)
p6 <- TSNEPlot(immune.combined, reduction = "tsne", split.by = "disease.ident",label = TRUE)
p7 <- DimPlot(immune.combined, reduction = "umap",label = TRUE)
p8 <- TSNEPlot(immune.combined, reduction = "tsne", label = TRUE)
grid.arrange(p7,p5,p8,p6, ncol = 2, nrow = 2)
p5|p6
# 
# # Find differentially expressed genes
# ## performing DE on the integrated assay is not statistically valid.
DefaultAssay(immune.combined) <- "RNA"
marker_0 <- FindConservedMarkers(immune.combined, 
                                 ident.1 = 0, 
                                 grouping.var = "disease.ident", 
                                 #assay="RNA",verbose = FALSE
                                 )
head(marker_0)

# let's visualize top features
FeaturePlot(immune.combined, features = c('CD22'), min.cutoff = 'q10')










# run Harmony -----------
cf2.harmony <- cf.2 %>%
  RunHarmony(group.by.vars = 'disease.ident', plot_convergence = FALSE)

cf2.harmony@reductions

cf2.harmony.embed <- Embeddings(cf2.harmony, "harmony")
cf2.harmony.embed[1:10,1:10]



# Do UMAP and clustering using ** Harmony embeddings instead of PCA **
cf2.harmony <- cf2.harmony %>%
  RunUMAP(reduction = 'harmony', dims = 1:20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20) %>%
  FindClusters(resolution = 1)

# visualize 
after <- DimPlot(cf2.harmony, reduction = 'umap', group.by = 'disease.ident',label = TRUE)
DimPlot(cf2.harmony, reduction = 'umap', label = TRUE)

DefaultAssay(cf2.harmony)

markers_h9 <- FindConservedMarkers(cf2.harmony,
                                         ident.1 = 9,
                                         grouping.var = 'disease.ident')

head(markers_h9)

# let's visualize top features
FeaturePlot(cf2.harmony, features = c('KRT19'), min.cutoff = 'q10')


# findMarkers between conditions ---------------------
cf2.harmony$celltype.cnd <- paste0(cf2.harmony$seurat_clusters,'_', cf2.harmony$disease.ident)
View(cf2.harmony@meta.data)
Idents(cf2.harmony) <- cf2.harmony$celltype.cnd

DimPlot(cf2.harmony, reduction = 'umap', label = TRUE)

# find markers
b.interferon.response <- FindMarkers(cf2.harmony, ident.1 = '9_CF', ident.2 = '9_CTRL')

head(b.interferon.response)

# plotting conserved features vs DE features between conditions
FeaturePlot(cf2.harmony, features = c('FABP4', 'AC105402.3', 'KRT19'), split.by = 'disease.ident', min.cutoff = 'q10')










saveRDS(immune.combined, file = "./results/cf2.combine_Jan18.rds")
saveRDS(cf.2, file = "./results/cf2_Jan18.rds")

















######mapping based on cell.type.ident
clusterdata_0 <- FetchData(cf.2, vars =  c("cell.type.ident", "seurat_clusters"))
clusterdata <- FetchData(immune.combined, vars =  c("cell.type.ident", "seurat_clusters"))
cluster_comp <-table(clusterdata_0$seurat_clusters,clusterdata$seurat_clusters)



tb_cluster <- table(clusterdata$cell.type.ident,clusterdata$seurat_clusters)
print(tb_cluster)
#based on the resutlts, we can define cluster 4-Mono,5-AM,7-Epithelial,8-(B,pDC,T_NK),0,1,2,3,6-(Multiplets,NG),cell type DC not sure


VlnPlot(immune.combined, features="CD22",split.by = 'disease.ident')
VlnPlot(immune.combined, features="CD22")

######
# tb_cluster_0 <- table(clusterdata_0$cell.type.ident,clusterdata_0$seurat_clusters)
# print(tb_cluster_0)
VlnPlot(cf.2, features="CD22")
######























################################################################################
################################################################################
#Subsetting CF.2 based on markers for Monocytes(Mono), Alveolar Macrophages(AM),
# & Monocyte derived Macrophages(MoMO)

#look for markers, remove clusters, redo downstream analysis, repeat

###############################################
##Subset 1

#CF.2 MARKERS dataframe used to manually label clusters based on highly
#expressed genes from literature

#common macrophage markers
VlnPlot(cf.2, features=c("FNIP2", "ALCAM"))

#remove neutrophil clusters and epitheleal cells
cf.2.small.1<-subset(cf.2, idents=c(9,7,5,4,3))   #3983 samples, 35451 features

#redo downstream with cf2.subset1
cf.2.small.1 <- FindVariableFeatures(cf.2.small.1, selection.method = "vst", nfeatures = 3000)

#scale and center data, regress out variables
cf.2.small.1<-ScaleData(cf.2.small.1, 
                        features=VariableFeatures(cf.2.small.1), verbose=F)
cf.2.small.1 <-ScaleData(cf.2.small.1, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)

#subset1 PCA dimensionaily reduction on subset of scaled variable features
cf.2.small.1 <- RunPCA(cf.2.small.1,
                       features = VariableFeatures(object = cf.2.small.1))

#visualize pc's
DimHeatmap(cf.2.small.1, dims = 1:20,cells=500, balanced=T, reduction = "pca")

#Determine "dimensionality" of subset1
ElbowPlot(cf.2.small.1, ndims=20, reduction="pca")

#######Cluster the cells
#nearest neighbors
cf.2.small.1 <- FindNeighbors(cf.2.small.1, dims = 1:12)

#cluster with Louvain (no resolution size outlined in methods)
cf.2.small.1 <- FindClusters(cf.2.small.1, resolution = 0.5)

##########Run non-linear dimensional reduction
#subset1 tsne
cf.2.small.1<-RunTSNE(cf.2.small.1, reduction="pca", dims=1:12)

#visualize
DimPlot(cf.2.small.1, reduction="tsne", label=T)

#find clusters with monocyte/macro upreg
VlnPlot(cf.2.small.1, features=c("SERPINB9", "TIMP1", "JARID2", #monocyte top
                                 "APOE","STMN1", "ANXA1",  #AM top
                                 "FNIP2", "ALCAM", "GPR183"      #M top
))

#compare cluster 2 expressions to 5 and 7
cluster2<-FindMarkers(cf.2.small.1, ident.1 =2, 
                      ident.2= c(5,7),logfc.threshold = 0.25,
                      test.use = "roc", only.pos = TRUE)

####################################################################
#SUBSET 2

#remove clusters with low expression of above 
cf.2.small.2<-subset(cf.2.small.1, idents=c(2,4,5,6,7
                                            #,8
                                            )) 
#35451 features, 1998 samples

#redo downstream with cf2.subset1
cf.2.small.2 <- FindVariableFeatures(cf.2.small.2, selection.method = "vst", nfeatures = 3000)

#scale to mean 0 and variance 1 & regress out percentage number of mitochondrial genes
cf.2.small.2<-ScaleData(cf.2.small.2, 
                        features=VariableFeatures(cf.2.small.2), verbose=F)
cf.2.small.2 <-ScaleData(cf.2.small.2, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)


#subset1 PCA dimensionaily reduction on subset 2 of scaled variable features
cf.2.small.2 <- RunPCA(cf.2.small.2, 
                       features = VariableFeatures(object = cf.2.small.2))

#Determine "dimensionality" of subset2
ElbowPlot(cf.2.small.2, ndims=20, reduction="pca")

#######Cluster the cells (subset2)
#nearest neighbors 
cf.2.small.2 <- FindNeighbors(cf.2.small.2, dims = 1:15) #higher dim number

#cluster with Louvain 
cf.2.small.2 <- FindClusters(cf.2.small.2, resolution = 0.5)

##########Run non-linear dimensional reduction
#subset2 tsne
cf.2.small.2<-RunTSNE(cf.2.small.2, reduction="pca", dims=1:15)

#visualize
DimPlot(cf.2.small.2, reduction="tsne", label=T)

#Get features from this subset
cf2.small.2.markers <- FindAllMarkers(cf.2.small.2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CF2.SMALL.2.MARKERS<-as.data.frame(cf2.small.2.markers %>%
                                     group_by(cluster) %>%
                                     slice_max(n = 10, order_by = avg_log2FC)
)

write.csv(CF2.SMALL.2.MARKERS, "CF_Data/small.2.markers(2).csv")
#########################################################################
##subset 3

#find clusters with monocyte/macro upreg
VlnPlot(cf.2.small.2, features=c("SERPINB9", "TIMP1", "JARID2", #monocyte top
                                 "APOE","STMN1", "ANXA1",  #AM top
                                 "FNIP2", "ALCAM", "GPR183"      #M top
))

#remove clusters based on top 10 markers(cf.2.small.2 df) and vln plot
cf.2.small.3<-subset(cf.2.small.2, idents=c(0,1,2,3,7))
#1698 samples

#redo downstream with cf2.subset3
cf.2.small.3 <- FindVariableFeatures(cf.2.small.3, selection.method = "vst", nfeatures = 3000)

#scale to mean 0 and variance 1 & regress out percentage number of mitochondrial genes
cf.2.small.3<-ScaleData(cf.2.small.3, features=VariableFeatures(cf.2.small.3), verbose=F)
cf.2.small.3 <-ScaleData(cf.2.small.3, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)

#subset1 PCA dimensionaily reduction on subset 3 of scaled variable features
cf.2.small.3 <- RunPCA(cf.2.small.3, 
                       features = VariableFeatures(object = cf.2.small.3))

#visualize pc dimensions
DimHeatmap(cf.2.small.3, dims=1:20, reduction="pca", cells=500,balanced = T)

#Determine "dimensionality" of subset1
ElbowPlot(cf.2.small.3, ndims=20, reduction="pca")

#######Cluster the cells (subset2)
#nearest neighbors 
cf.2.small.3 <- FindNeighbors(cf.2.small.3, dims = 1:15) #higher dimension number

#cluster with Louvain
cf.2.small.3 <- FindClusters(cf.2.small.3, resolution = 0.5)

##########Run non-linear dimensional reduction(UMAP/tSNE)
#subset1 tsne
cf.2.small.3<-RunTSNE(cf.2.small.3, reduction="pca", dims=1:15)

#visualize
DimPlot(cf.2.small.3, reduction="tsne", label=T)

#Get features from this subset
cf2.small.3.markers <- FindAllMarkers(cf.2.small.3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CF2.SMALL.3.MARKERS<-as.data.frame(cf2.small.3.markers %>%
                                     group_by(cluster) %>%
                                     slice_max(n = 10, order_by = avg_log2FC)
)

write.csv(CF2.SMALL.3.MARKERS, "CF_Data/cf.2.small.3.markers.csv")

small.3.topgenes<-cf2.small.3.markers %>%
  top_n(n = 10, wt = -(p_val_adj)) 

########################
#Rename clusters for subset 3(based on literature & top 10 genes)
new.cluster.ids <- c("Mono1", "Mono2", "AM1", "Mono3","AM2", "Mult")
names(new.cluster.ids) <- levels(cf.2.small.3)
cf.2.small.3 <- RenameIdents(cf.2.small.3, new.cluster.ids)

#replot tsne with new names
DimPlot(cf.2.small.3, reduction="tsne", label=T)

#compare cf vs. control subset 3
VlnPlot(cf.2.small.3, 
        features = small.3.topgenes$gene , 
        split.by = "disease.ident",split.plot=T,log = T)+
  theme(legend.position = "right")

###########################################################
#look at known necroptosis and senescence markers in subset 3
#necrop: RIP1, RIP2, MLKL
#senesc: CDKN2A, CDKN1A,TP53

VlnPlot(cf.2.small.3, 
        features=c("RIPK1", "RIPK2", "MLKL", "CDKN2A", "CDKN1A", "TP53"), 
        log=T)
#three Monocyte clusters have markers present

#compare ripk2 and cdkn1a cf/control that are present
VlnPlot(cf.2.small.3, features=c("RIPK2", "CDKN1A"), split.by="disease.ident",
        split.plot = T, log=T)+
  theme(legend.position = "right")

#############################################################
#subset 4
#remove cluster with multiple cell types

cf.2.small.4<-subset(cf.2.small.3, idents=c("Mono1", "Mono2", "Mono3", 
                                            "AM1","AM2")) #1634 samples

#redo downstream with cf2.subset4
cf.2.small.4<- FindVariableFeatures(cf.2.small.4, selection.method = "vst", nfeatures = 3000)

#scale to mean 0 and variance 1 & regress out percentage number of mitochondrial genes

cf.2.small.4<-ScaleData(cf.2.small.4, 
                        features=VariableFeatures(cf.2.small.4), verbose=F)
cf.2.small.4 <-ScaleData(cf.2.small.4, 
                         vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), 
                         do.scale=F, do.center=F, verbose=F)


#subset4 PCA dimensionaily reduction on subset 4 of scaled variable features
cf.2.small.4 <- RunPCA(cf.2.small.4, features = VariableFeatures(object = cf.2.small.4))

DimHeatmap(cf.2.small.4, dims=1:20, cells=500, balanced=T, reduction = "pca")

#Determine "dimensionality" of subset4
ElbowPlot(cf.2.small.4, ndims=20, reduction="pca")

#######Cluster the cells (subset4)
#nearest neighbors 
cf.2.small.4 <- FindNeighbors(cf.2.small.4, dims = 1:12) 

#cluster with Louvain
cf.2.small.4 <- FindClusters(cf.2.small.4, resolution = 0.8)

##########Run non-linear dimensional reduction(UMAP/tSNE)
#subset4 tsne
cf.2.small.4<-RunTSNE(cf.2.small.4, reduction="pca", dims=1:12)

#visualize
DimPlot(cf.2.small.4, reduction="tsne", label=T)

DimPlot(cf.2.small.4, reduction="pca", label=T)

DimHeatmap(cf.2.small.4, reduction = "pca", dims = 1:12,cells = 500, balanced=T)

#Get features from this subset
cf2.small.4.markers <- FindAllMarkers(cf.2.small.4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
CF2.SMALL.4.MARKERS<-as.data.frame(cf2.small.4.markers %>%
                                     group_by(cluster) %>%
                                     slice_max(n = 20, order_by = avg_log2FC)
)

write.csv(CF2.SMALL.4.MARKERS, "CF_Data/cf.2.small.4.markers.csv")

small.4.topgenes<-cf2.small.4.markers %>%
  top_n(n = 10, wt = -(p_val_adj)) 

#compare cf vs. control for top 10 genes
VlnPlot(cf.2.small.4, 
        features = small.4.topgenes$gene , 
        split.by = "disease.ident",split.plot=T,log = T)+
  theme(legend.position = "right")

########################
#Rename clusters for subset 4(literature and top 20 DF)
new.cluster.ids <- c("Mono1", "AM1", "Mono2/MoMO", "Mono3","Mono4", "AM2")
names(new.cluster.ids) <- levels(cf.2.small.4)
cf.2.small.4 <- RenameIdents(cf.2.small.4, new.cluster.ids)

#replot tsne with new names
DimPlot(cf.2.small.4, reduction="tsne",shape.by = "disease.ident",label=T)


#look at known necroptosis and senescence markers in subset 4
#necrop: RIP1, RIP2, MLKL
#senesc: CDKN2A, CDKN1A,TP53

VlnPlot(cf.2.small.4, 
        features=c("RIPK1", "RIPK2", "MLKL", "CDKN2A", "CDKN1A", "TP53"), 
        log=T)
#all monocytes have expression upreg

#compare ripk2 and cdkn1a cf/control that are present
VlnPlot(cf.2.small.4, features=c("RIPK2", "CDKN1A"), split.by="disease.ident",
        split.plot = T, log=T)+
  theme(legend.position = "right")

saveRDS(cf.2.small.4, file="CF_Data/cf.2.small.4.rds")

######################################################
#try different downstream dim=8 (no change)
{
  #nearest neighbors 
  cf.2.small.4.8 <- FindNeighbors(cf.2.small.4, dims = 1:8) 
  
  #cluster with Louvain
  cf.2.small.4.8 <- FindClusters(cf.2.small.4.8, resolution = 0.8)
  
  #subset4 tsne
  cf.2.small.4.8<-RunTSNE(cf.2.small.4.8, reduction="pca", dims=1:8)
  #visualize
  DimPlot(cf.2.small.4.8, reduction="tsne", label=T)
  
  DimHeatmap(cf.2.small.4.8, reduction = "pca", dims = 1:8,
             cells = 500, balanced=T)
}
#######################################################
###Markers from Yong
#no MHCII, cd64, CD11b, ly6c, CD169
#look at all cf.2
VlnPlot(cf.2, 
        features=c("RIPK3" ,"ADGRE1", "CD68", 
                   "MERTK", "CD14"), log=T)

#look at subset
VlnPlot(cf.2.small.4, 
        features=c("RIPK3" ,"ADGRE1", "CD68", 
                   "MERTK", "CD14"), log=T)















