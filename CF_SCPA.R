library(SCPA)
library(msigdbr)
library(Seurat)
library(dplyr)
library(ggplot2)


# look at the gene sets in GSEA
# print(msigdbr_collections(),n=23)

# read and plot rds-------------------------
CF_integration.1 <- readRDS("./results/CF_integration_1.rds")
Idents(CF_integration.1) <- CF_integration.1@meta.data$rpca_clusters_0.1
DimPlot(CF_integration.1)

# get the cluster subest --------------------------
# # combine cluster 1 and cluster 5 (AM)------------------
# scpa_seurat_1 <- subset(CF_integration.1, idents = 1)
# scpa_seurat_5 <- subset(CF_integration.1, idents = 5)
# cf_scpa_1_5 <- merge(x = scpa_seurat_1, y = scpa_seurat_5)
# cf_scpa_1_5 <- JoinLayers(cf_scpa_1_5)
# scpa_data_1_5 <- GetAssayData(cf_scpa_1_5, assay = "RNA", layer = "data")
# scpa_data_1_5 <- as.matrix(scpa_data_1_5)

# a better way
cf_scpa_1_5 <- subset(CF_integration.1, idents = c(1,5), invert = FALSE)
scpa_data_1_5 <- GetAssayData(cf_scpa_1_5, assay = "RNA", layer = "data")
scpa_data_1_5 <- as.matrix(scpa_data_1_5)
# cluster 2---------------------------------------------
scpa_seurat_2 <- subset(CF_integration.1, idents = 2)
scpa_data_2 <- GetAssayData(scpa_seurat_2, assay = "RNA", layer = "data")
scpa_data_2 <- as.matrix(scpa_data_2)

# scpa package can also compare seurat object, but no sure how to do this for combining 2 object

# get the pathway gene set-----------------------------------------
### homo sapiens is human
pathways <- msigdbr(species ="Homo sapiens",category ="C5"
                    ,subcategory ="GO:BP") %>%
  format_pathways()
head(pathways)

# pathway analysis----------------------------------------------------
scpa_out <- compare_pathways(samples = list(scpa_data_1_5,scpa_data_2), 
                             pathways = pathways)
#scpa_data_2:MDM --- CF
head(scpa_out,20)

write.csv(scpa_out,file="./results/SCPA/SCPAout.csv")
# a higher qval translates to larger pathway differences between conditions
# a fold change (FC) enrichment score will also be calculated. 
#Instead of looking at pathway enrichment, 
#SCPA assesses changes in the multivariate distribution of a pathway. 
#However, pathways that show enrichment in a given population will also necessarily show large changes in multivariate distribution. 
#This means that with SCPA you’re able to detect 1) enriched pathways and 2) non-enriched pathways that have transcriptional changes that are independent of enrichment. 
#This output is reflected in the qval, and we recommend people to use this as their primary statistic. 
#Here, the larger the qval, the larger the change in pathway ‘activity’. 
#As SCPA measures multivariate distributions, there will be pathways 
#that show significantly large qvals, but no overall fold change/enrichment in a given population. 
#Whilst these pathways are not enriched, we know that these distribution changes identified by SCPA are still important for cellular behaviour.


scpa_out2 <- compare_pathways(samples = list(scpa_data_2,scpa_data_1_5), 
                             pathways = pathways)
#scpa_data_2:MDM --- CF
head(scpa_out2,20)
write.csv(scpa_out,file="./results/SCPA/SCPAout2.csv")


# plot------------
plot_rank(scpa_out = scpa_out, 
          pathway = "lymphocyte", 
          base_point_size = 2, 
          highlight_point_size = 4)

plot_heatmap(scpa_out, 
             highlight_pathways = "lymphocyte",
             show_row_names = F)

scpa_out <- scpa_out %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

aa_path <- scpa_out %>% 
  filter(grepl(pattern = "reactome_arachi", ignore.case = T, x = Pathway))

ggplot(scpa_out, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = scpa_out$color, stroke = 0.3) +
  geom_point(data = aa_path, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)




scpa_out2 <- scpa_out2 %>%
  mutate(color = case_when(FC > 5 & adjPval < 0.01 ~ '#6dbf88',
                           FC < 5 & FC > -5 & adjPval < 0.01 ~ '#84b0f0',
                           FC < -5 & adjPval < 0.01 ~ 'mediumseagreen',
                           FC < 5 & FC > -5 & adjPval > 0.01 ~ 'black'))

aa_path2 <- scpa_out2 %>% 
  filter(grepl(pattern = "reactome_arachi", ignore.case = T, x = Pathway))

ggplot(scpa_out2, aes(-FC, qval)) +
  geom_vline(xintercept = c(-5, 5), linetype = "dashed", col = 'black', lwd = 0.3) +
  geom_point(cex = 2.6, shape = 21, fill = scpa_out2$color, stroke = 0.3) +
  geom_point(data = aa_path2, shape = 21, cex = 2.8, fill = "orangered2", color = "black", stroke = 0.3) +
  xlim(-20, 80) +
  ylim(0, 11) +
  xlab("Enrichment") +
  ylab("Qval") +
  theme(panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        aspect.ratio = 1)




plot <- plot_heatmap(scpa_out2, 
                 highlight_pathways = "LYMPHOCYTE",
                 column_names = "MDM vs AM",
                 show_row_names = F)




# save plots in PDF--------------------
#open PDF
pdf(file="./results/SCPA/plots_scpa.pdf")

#save plots to PDF
plot

#turn off PDF plotting
dev.off() 
