setwd("/Users/menglinghe/Library/Mobile Documents/com~apple~CloudDocs/UTK/GRA-UTK/UTK_R/CF_project")

# load libraries
library(dplyr)
library(tidyverse)
library(Seurat)
library(patchwork)
library(data.table)
library(R.utils)
library(ggplot2)
library(patchwork)
library(gridExtra)
library(harmony)
set.seed(412)


# read and plot rds-------------------------
cf.1 <- readRDS("./results/cf_1.rds") # will use this rds for further analysis
CF_integration.1 <- readRDS("./results/CF_integration_1.rds") # will use this clustering to do subset
# use integration method=rpca, resolution=0.1
Idents(cf.1) <- CF_integration.1@meta.data$rpca_clusters_0.1
DimPlot(cf.1)


# subset of cluster 0,2,1,5----------------------------------
subset_cf1 <- subset(cf.1, idents = c(0,1,2,5), invert = FALSE)
DimPlot(subset_cf1)




# data processing: FindVariableFeatures and scale----------------------------
# feature selection(find highly variable features)
subset_cf1 <- FindVariableFeatures(subset_cf1, selection.method = "vst", nfeatures = 3000)
# scale to mean 0 and variance 1
all.genes <- rownames(subset_cf1)
subset_cf1<-ScaleData(subset_cf1, features=all.genes, verbose=F)
# regress out percentage number of mitochondrial genes
subset_cf1 <-ScaleData(subset_cf1, vars.to.regress = c("nFeature_RNA", "nCount_RNA","percent.mt"), do.scale=F, do.center=F, verbose=F)



# Perform streamlined (one-line) integrative analysis--------------------
CFsubset_integration <- IntegrateLayers(
  object = subset_cf1, method = RPCAIntegration,
  orig.reduction = "pca", new.reduction = 'integrated.rpca',
  verbose = FALSE)

CFsubset_integration <- FindNeighbors(CFsubset_integration, reduction = "integrated.rpca", dims = 1:15)
CFsubset_integration <- FindClusters(CFsubset_integration, resolution =c(0.1,0.3,0.5), cluster.name = c("subcluster_rpca0.1","subcluster_rpca0.3","subcluster_rpca0.5"))
CFsubset_integration <- FindClusters(CFsubset_integration, resolution =c(0.7,1,1.2), cluster.name = c("subcluster_rpca0.7","subcluster_rpca1","subcluster_rpca1.2"))
CFsubset_integration <- RunUMAP(CFsubset_integration, reduction = "integrated.rpca", dims = 1:15, reduction.name = "umap.rpca")
p2 <- DimPlot(
  CFsubset_integration,
  reduction = "umap.rpca",
  group.by = c("disease.ident", "cell.type.ident"),
  combine = FALSE, label.size = 2,label = T
)


#  find the resolution ------------------------
DimPlot(CFsubset_integration, reduction = "umap.rpca"
        , group.by = c("subcluster_rpca0.1","subcluster_rpca0.3","subcluster_rpca0.5","subcluster_rpca0.7","subcluster_rpca1","subcluster_rpca1.2")
        ,combine = TRUE
        ,label = T)
# the clustering is similar to non subset clustering when resolution=0.1
# 
# DimPlot(CFsubset_integration, reduction = "umap.unintegrated",label = T)|DimPlot(
#   CFsubset_integration, reduction = "umap.unintegrated", group.by = "subcluster_rpca0.3",label = T)




# choose resolution = 0.1 and visualization
Idents(CFsubset_integration) <- Idents(subset_cf1)
# compare the subclustering with the original clustering
DimPlot(CFsubset_integration, reduction = "umap.rpca",label = T)|DimPlot(
  CFsubset_integration, reduction = "umap.rpca", group.by = "subcluster_rpca0.1",label = T)

# visulize the subset clustering in CF and CTRL group
p1 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,group.by = "subcluster_rpca0.1",label = T)
p2 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,group.by = "subcluster_rpca0.1",split.by = "disease.ident",label = T)
grid.arrange(p1,p2,nrow=2)





# choose resolution = 0.5 and visualization
Idents(CFsubset_integration) <- Idents(subset_cf1)
# compare the subclustering with the original clustering
DimPlot(CFsubset_integration, reduction = "umap.rpca",label = T)|DimPlot(
  CFsubset_integration, reduction = "umap.rpca", group.by = "subcluster_rpca0.5",label = T)

# visulize the subset clustering in CF and CTRL group
p1 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,group.by = "subcluster_rpca0.5",label = T)
p2 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,group.by = "subcluster_rpca0.5",split.by = "disease.ident",label = T)
grid.arrange(p1,p2,nrow=2)



# choose resolution = 0.7 and visualization
Idents(CFsubset_integration) <- Idents(subset_cf1)
# compare the subclustering with the original clustering
DimPlot(CFsubset_integration, reduction = "umap.rpca",label = T)|DimPlot(
  CFsubset_integration, reduction = "umap.rpca", group.by = "subcluster_rpca0.7",label = T)

# visulize the subset clustering in CF and CTRL group
p1 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,group.by = "subcluster_rpca0.7",label = T)
p2 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
        ,group.by = "subcluster_rpca0.7",split.by = "disease.ident",label = T)
grid.arrange(p1,p2,nrow=2)





# rename the subgroup (rpca0.1)------------------------------------
Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.1
newnames0.1 <- c('PMN','AMs_1','MDM','AMs_2')
names(newnames0.1) <- as.character(sort(as.numeric(levels(CFsubset_integration))))
CFsubset_integration <- RenameIdents(CFsubset_integration,newnames0.1)

# visulize the subset clustering in CF and CTRL group in renamed clustering
p1 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,label = T,label.size = 3)+ NoLegend()
p2 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,split.by = "disease.ident",label = T,label.size = 3)+ NoLegend()
grid.arrange(p1,p2,nrow=2)




# rename the subgroup (rpca0.5)------------------------------------
Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.5
# g0 <- c('1','2','4','8','9')#PMN
# g1 <- as.character(c(0,5,7,10))#AM
# g2 <-as.character(c(3,6,11))#MDM
newnames0.5 <- c('AMs_1','PMNs_1','PMNs_2','MDMs_1','PMNs_3','AMs_2','MDMs_2','AMs_3','PMNs_4','PMNs_5'
              ,'AMs_4','MDMs_2')
names(newnames0.5) <- as.character(sort(as.numeric(levels(CFsubset_integration))))
CFsubset_integration <- RenameIdents(CFsubset_integration,newnames0.5)

# visulize the subset clustering in CF and CTRL group in renamed clustering
p1 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,label = T,label.size = 3)+ NoLegend()
p2 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
              ,split.by = "disease.ident",label = T,label.size = 3)+ NoLegend()
grid.arrange(p1,p2,nrow=2)





# 
# 
# # rename the subgroup (rpca0.7)------------------------------------
# Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.7
# # g0 <- c(0,'1','2','8','10','12')#PMN
# # g1 <- c('4','5','6','7','13')#AM
# # g2 <-as.character(c(3,9,11,14))#MDM
# newnames0.7 <- c('PMNs_1','PMNs_2','PMNs_3','MDMs_1','AMs_1','AMs_2'
#               ,'AMs_3','AMs_4','PMNs_4','MDMs_2','PMNs_5'
#               ,'MDMs_3','PMNs_6','AMs_5','MDMs_4')
# names(newnames0.7) <- as.character(sort(as.numeric(levels(CFsubset_integration))))
# CFsubset_integration <- RenameIdents(CFsubset_integration,newnames0.7)
# 
# # visulize the subset clustering in CF and CTRL group in renamed clustering
# p1 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
#               ,label = T,label.size = 3)+ NoLegend()
# p2 <- DimPlot(CFsubset_integration, reduction = "umap.rpca"
#               ,split.by = "disease.ident",label = T,label.size = 3)+ NoLegend()
# grid.arrange(p1,p2,nrow=2)






# number of samples in each cluster----------------------------------------
CFsub_summary2 <- table( CFsubset_integration@meta.data$disease.ident,Idents(CFsubset_integration))
tab <- cbind(CFsub_summary2, Total = rowSums(CFsub_summary2))
print(tab)


# Barplot of proportion of cells in each cluster by sample
ggplot(CFsubset_integration@meta.data) +
  geom_bar(aes(x=subcluster_rpca0.5, fill=disease.ident), position=position_fill()) 

# 
# # cell by cluster ------------------------
# cluster_bycell_subset <- CFsubset_integration@meta.data["subcluster_rpca0.3"]
# cluster_bycell <- read.csv("./results/cluster_bycell.csv",row.names=1)
# cluster_bycell_subset0 <- cluster_bycell %>% filter(cluster_bycell$rpca_clusters_0.1 %in% c('0', '1', '2','5'))
# cluster_bycell_compare <- merge(cluster_bycell_subset0, cluster_bycell_subset, 
#                           by = 'row.names', all = TRUE)
# 
# write.csv(cluster_bycell_compare, "./results/cluster_bycell_compare.csv", row.names=TRUE)










# Pre finding markers------------------------
#CFsubset_integration <- JoinLayers(CFsubset_integration)
#saveRDS(CFsubset_integration, file = "./results/CFsubset_integration.rds")
Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.5
n_cluster <- length(levels(Idents(CFsubset_integration)))









# # Find all markers------------------------------------------
# markers_findall <- FindAllMarkers(CFsubset_integration, only.pos = TRUE, min.pct = 0.25,
#                               logfc.threshold = 0.25)
# markers_findall<-as.data.frame(markers_findall %>%
#                              group_by(cluster) %>%
#                              slice_max(n = 20, order_by = avg_log2FC)
# )
# 
# 
# 
# 
# 
# marker_top5 <- markers_findall%>%group_by(cluster) %>% top_n(n=5,wt = avg_log2FC)
# 
# options(repr.plot.width = 13, repr.plot.height=13)
# DoHeatmap(CFsubset_integration, features = marker_top5$gene, disp.min=-2.5, disp.max=2.5)
# 
# 
# 
# # visualize MDMs_2 (cluster 9) just an example
# features2 <- markers_findall[markers_findall$cluster=="MDMs_2",]$gene
# FeaturePlot(CFsubset_integration, features = features2[1:3], split.by = 'disease.ident', min.cutoff = 'q10',reduction = 'umap.rpca')
# CFsubset_integration <- RenameIdents(CFsubset_integration,newnames)
# DoHeatmap(CFsubset_integration, features = features2[1:3], size = 3, angle = 0,  draw.lines = T)
# 
# 







# find conserve markers (see 'subsetanalysis_cf1_resolution' for a clean version)--------------------------------------
# find marker
marker = list()
for (i in 1:n_cluster){
  marker[[i]] <- FindConservedMarkers(CFsubset_integration, ident.1 = i-1
                                      , grouping.var = "disease.ident"
                                      #, only.pos = TRUE
                                      ,logfc.threshold = 0.25
                                      ,verbose = FALSE)
}
#head(marker[[13]],2)

# rank and get the first 20 markers (firstly filer p value less than 5% and than rank by fold change)
N=20
top_marker <- list()

for (i in c(1:9,11:12)){
  top_marker[[i]] <- marker[[i]] %>% filter(max_pval < 0.05) %>%
    mutate(avg_fc = (CF_avg_log2FC + CTRL_avg_log2FC) /2) %>% 
    slice_max(n = N, order_by = avg_fc)
    # top_n(n = N, 
    #       wt = avg_fc),
}

# for cluster including only one group
for (i in 10){
  top_marker[[i]] <- marker[[i]] %>% filter(.[[5]] < 0.05) %>%
    slice_max(n = N, order_by = CF_avg_log2FC)
}


# convert marker list into df
top_marker <- lapply(top_marker, function(df) rownames_to_column(df,var = "gene"))

marker.df <- bind_rows(top_marker, .id = 'cluster')
marker.df$cluster <- as.character(as.integer(marker.df$cluster)-1)
write.csv(marker.df, "./results/subcluster_marker0.5.csv", row.names=TRUE)







# plot the heatmap of top 10 markers
consmarker_top10 <- marker.df%>%group_by(cluster) %>% top_n(n=10)
CFsubset_integration <- RenameIdents(CFsubset_integration,newnames)
options(repr.plot.width = 13, repr.plot.height=13)
hp <- DoHeatmap(CFsubset_integration, features = consmarker_top10$gene
          ,size =3.5
          , disp.min=-2.5, disp.max=2.5)
hp+ ggtitle("Gene marker heatmap")+theme(plot.title = element_text(color="red", size=12))



hp_cf <- DoHeatmap(subset(x = CFsubset_integration,subset = disease.ident =="CF")
          , features = consmarker_top10$gene
          , size =3.5
          , disp.min=-2.5, disp.max=2.5)
hp_cf+ ggtitle("Gene marker heatmap -- CF")+theme(plot.title = element_text(color="red", size=12))


hp_ctrl<- DoHeatmap(subset(x = CFsubset_integration,subset = disease.ident =="CTRL")
          , features = consmarker_top10$gene
          , size =3.5
          , disp.min=-2.5, disp.max=2.5) 
hp_ctrl+ ggtitle("Gene marker heatmap -- CTRL")+theme(plot.title = element_text(color="red", size=12))

# 
# png('./results/heatmap1.png')
# hp+ ggtitle("Gene marker heatmap")+theme(plot.title = element_text(color="red", size=12))
# png('./results/heatmap_cf.png')
# hp_cf+ ggtitle("Gene marker heatmap -- CF")+theme(plot.title = element_text(color="red", size=12))
# png('./results/heatmap_ctrl.png')
# hp_ctrl+ ggtitle("Gene marker heatmap -- CTRL")+theme(plot.title = element_text(color="red", size=12))
# dev.off()






# 
# # an analysis of the top genes considering the ranking method-------------------------------------
# View(top_marker[[1]])
# View(top_marker[[2]])
# FeaturePlot(CFsubset_integration, features=head(top_marker[[1]],4)$gene, reduction ="umap.rpca", min.cutoff = 'q10')
# VlnPlot(CFsubset_integration, features=head(top_marker[[1]],6)$gene)
# # for resolution=0.7, there are 15 clusters, and choose the top markers by sorting the average fold change
# # the feature plot and vln plot is not very consistent
# # this is because in each cluster, the number of cells from CF and CTRL is very unbalanced, 
# # so su=imply take the average is not very meaningful
# # for example, here in cluster 0, majority of the cells are from CF. and for gene "PADI4" the fold change of CF group is 0.66 which is very small.
# # but the the fold change of Ctrl group is very large 7.31 which makes it the first gene marker in cluster 0
# # however, since ctrl is only a small number, so the expression is not very high from vln plot and feature plot
# FeaturePlot(CFsubset_integration, features=rownames(head(top_wmarker0,4)), reduction ="umap.rpca", min.cutoff = 'q10')
# VlnPlot(CFsubset_integration, features=rownames(head(top_wmarker0,4)))
# 
# 











# findMarkers between conditions ---------------------
n_cluster <- length(levels(CFsubset_integration$subcluster_rpca0.5))
CFsubset_integration$cluster.cnd <- paste0(CFsubset_integration$subcluster_rpca0.5,'.', CFsubset_integration$disease.ident)
print(length(levels(Idents(CFsubset_integration))))
print(n_cluster)
Idents(CFsubset_integration) <- CFsubset_integration$cluster.cnd
DimPlot(CFsubset_integration, reduction = 'umap.rpca', label = TRUE)


# find marker
interferon_list = list()
for (i in c(1:9,11:12)){
  ident1 <- paste0(i-1,'.','CF')
  ident2 <- paste0(i-1,'.','CTRL')
  interferon_list[[i]] <- FindMarkers(CFsubset_integration, ident.1 = ident1, ident.2 = ident2)
}
head(interferon_list[[1]],5)

# convert marker list into df
interferon_list_top = list()
for (i in c(1:9,11:12)){
  interferon_list_top[[i]] <- interferon_list[[i]] %>% filter(p_val_adj < 0.05) %>%
    slice_max(n = N, order_by = avg_log2FC)
}
interferon_list_top[[10]] <- data.frame()


# convert marker list into df
interferon_list_top <- lapply(interferon_list_top, function(df) rownames_to_column(df,var = "gene"))

interferon_top <- bind_rows(interferon_list_top, .id = 'cluster')
interferon_top$cluster <- as.character(as.integer(interferon_top$cluster)-1)

write.csv(interferon_top, "./results/withinsubcluster_marker_0.5.csv", row.names=TRUE)


# plot the heatmap of top 10 markers
intermarker_top10 <- interferon_top%>%group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
# Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.5
# CFsubset_integration <- RenameIdents(CFsubset_integration,newnames)
Idents(CFsubset_integration) <- CFsubset_integration@meta.data$disease.ident
options(repr.plot.width = 13, repr.plot.height=13)
hp <- DoHeatmap(CFsubset_integration, features = intermarker_top10$gene
                ,size =3.5
                , disp.min=-2.5, disp.max=2.5)
hp+ ggtitle("Gene marker heatmap")+theme(plot.title = element_text(color="red", size=12))

FeaturePlot(CFsubset_integration, features = c('G0S2', 'PRKCB', 'ANO5','AC026369.3'), split.by = 'disease.ident', min.cutoff = 'q10')


