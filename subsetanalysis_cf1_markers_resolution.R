# Create function to get conserved markers for any given cluster------------
get_conserved <- function(cluster){
  FindConservedMarkers(CFsubset_integration,
                       ident.1 = cluster,
                       grouping.var = "disease.ident",
                       only.pos = TRUE) %>%
    rownames_to_column(var = "gene") %>%
    cbind(cluster_id = cluster, .)
}













# resolution =0.1-------------------
Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.1
n_cluster <- length(levels(Idents(CFsubset_integration)))
# 
CFsub_summary0.1 <- table( CFsubset_integration@meta.data$disease.ident,Idents(CFsubset_integration))
tab_0.1 <- rbind(CFsub_summary0.1, Total = colSums(CFsub_summary0.1))
print(tab_0.1)

# Iterate function across desired clusters
conserved_markers0.1 <- map_dfr(0:(n_cluster-1), get_conserved)

# Extract top 10 markers per cluster(try 3 different ranking)
# only filter p value for those clusters with cells from 2 samples
top10_0.1 <- conserved_markers0.1 %>% 
  filter(max_pval < 0.05)%>% 
  mutate(c =as.numeric(cluster_id+1)) %>% 
  mutate(avg_fc = (CF_avg_log2FC*tab_0.1[1,c] + CTRL_avg_log2FC*tab_0.1[2,c]) 
         /tab_0.1[3,c]) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

top10_0.1_2 <- conserved_markers0.1 %>% 
  filter(max_pval < 0.05)%>% 
  mutate(c =as.numeric(cluster_id+1)) %>% 
  mutate(avg_fc = (CF_avg_log2FC+ CTRL_avg_log2FC)/2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)

top10_0.1_3 <- conserved_markers0.1 %>% 
  filter(max_pval < 0.05)%>% 
  mutate(c =as.numeric(cluster_id+1)) %>% 
  mutate(avg_fc = log(((2^CF_avg_log2FC)*tab_0.1[1,c] + (2^CTRL_avg_log2FC)*tab_0.1[2,c]) 
         /tab_0.1[3,c]),2) %>% 
  group_by(cluster_id) %>% 
  top_n(n = 10, 
        wt = avg_fc)



DimPlot(CFsubset_integration, reduction = "umap.rpca",label = T)
FeaturePlot(CFsubset_integration, features=head(top10_0.1,4)$gene, reduction ="umap.rpca", min.cutoff = 'q10')
FeaturePlot(CFsubset_integration, features=head(top10_0.1_2,4)$gene, reduction ="umap.rpca", min.cutoff = 'q10')
FeaturePlot(CFsubset_integration, features=head(top10_0.1_3,4)$gene, reduction ="umap.rpca", min.cutoff = 'q10')

VlnPlot(CFsubset_integration, features=head(top10_0.1,6)$gene)
VlnPlot(CFsubset_integration, features=head(top10_0.1_2,6)$gene)
VlnPlot(CFsubset_integration, features=head(top10_0.1_3,6)$gene)

head(top10_0.1,10)$gene
head(top10_0.1_3,10)$gene
head(top10_0.1_2,10)$gene

options(repr.plot.width = 13, repr.plot.height=13)
hp <- DoHeatmap(CFsubset_integration
                , features = top10_0.1$gene
                ,size =3.5
                , disp.min=-2.5, disp.max=2.5)
hp+ ggtitle("Gene marker heatmap")+theme(plot.title = element_text(color="red", size=12))

DoHeatmap(CFsubset_integration
          , features = top10_0.1_2$gene
          ,size =3.5
          , disp.min=-2.5, disp.max=2.5)

DoHeatmap(CFsubset_integration
          , features = top10_0.1_3$gene
          ,size =3.5
          , disp.min=-2.5, disp.max=2.5)



# resolution =0.5-------------------
Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.5
n_cluster <- length(levels(Idents(CFsubset_integration)))
# 
CFsub_summary0.5 <- table( CFsubset_integration@meta.data$disease.ident,Idents(CFsubset_integration))
tab_0.5 <- rbind(CFsub_summary0.5, Total = colSums(CFsub_summary0.5))

# Iterate function across desired clusters
conserved_markers0.5 <- map_dfr(0:(n_cluster-1), get_conserved)



d <-conserved_markers0.5[conserved_markers0.5$cluster_id=='9',]


# resolution =0.7-------------------
Idents(CFsubset_integration) <- CFsubset_integration@meta.data$subcluster_rpca0.7
n_cluster <- length(levels(Idents(CFsubset_integration)))
# 
CFsub_summary0.7 <- table( CFsubset_integration@meta.data$disease.ident,Idents(CFsubset_integration))
tab_0.7 <- rbind(CFsub_summary0.7, Total = colSums(CFsub_summary0.7))

# Iterate function across desired clusters
conserved_markers0.7 <- map_dfr(0:(n_cluster-1), get_conserved)


write.csv(conserved_markers0.7, "./results/subclustering/subset_rpca0.7/conserved_markers0.7.csv", row.names=FALSE)
write.csv(conserved_markers0.5, "./results/subclustering/subset_rpca0.5/conserved_markers0.5.csv", row.names=FALSE)
write.csv(conserved_markers0.1, "./results/subclustering/subset_rpca0.1/conserved_markers0.1.csv", row.names=FALSE)


