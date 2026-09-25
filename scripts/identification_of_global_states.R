library(Seurat)
library(dplyr)

#Set Parameters
input.dir <- '/path_to_data/'
output.dir <- '/path_to_data/'
letters <- c("A", "B", "C", "D", "E", "G", "I", "K", "M")
samples <- c("sampleA", "sampleB", "sampleC", "sampleD", "sampleE", "sampleG", "sampleI", 
             "sampleK", "sampleM")
size <- c(rep(10500, 4), 5500, 3000, 4000, 5500, 3500)
cols <- c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Interface' = 'grey', 
          "Active" = '#FF3C31', "Passive" = "#5B8AFD", 'Intermediate' = '#b5179e')


##### 1. Process object and add metadata #########################################
#load object
load("/path_to_data/bc.hallmark.pathways.RData")

#plot object 
SpatialFeaturePlot(bc_hallmark_pathways, "NOTCH2", images = "sampleA", pt.size.factor  = size[1], interactive = F)

#remove NATs from object 
object <- UpdateSeuratObject(bc_hallmark_pathways)
object <- subset(object, subset = !(orig.ident %in% c('sample_F','sample_H', 'sample_J', 'sample_L', 'sample_N')))
SpatialFeaturePlot(object, "NOTCH2", images = "sampleA", pt.size.factor  = size[1])
rm(bc_hallmark_pathways)

#Add metadata
meta.data <- c()
for(sample in letters){
  border.anno <- read.csv(paste0('/path_to_data/Annotations/', 
                                 sample, '_border_annotations.csv'))
  meta.data <- rbind(meta.data, border.anno[,c("cell_types", 'final_annotations_w_activity',"final_annotations", "interface_w_border")])
}
object <- AddMetaData(object, meta.data)

#Check metadata
sample.num <- 8
SpatialPlot(object, 'final_annotations_w_activity',images = samples[sample.num], pt.size.factor  = size[sample.num], cols= cols)

#Visualize umap
object@meta.data$border.zone <- object$interface_w_border
object@meta.data[object$interface_w_border == 'Border',]$interface.zone <- object@meta.data[object$interface_w_border == 'Border',]$final_annotations_w_activity
DimPlot(object , group.by = 'final_annotations_w_activity', cols= cols)


#### 2. Extract and process Interface Object ########################################################

#Subset
object.inter <- subset(object, subset = final_annotations_w_activity %in% c('High', 'Medium', 'Low'))

#Get Umap
DefaultAssay(object.inter) <- "integrated"
object.inter <- FindVariableFeatures(object.inter)
object.inter <- RunPCA(object.inter, verbose = FALSE)
object.inter <- FindNeighbors(object.inter, reduction = "pca", dims = 1:30)
object.inter <- FindClusters(object.inter, verbose = FALSE)
object.inter <- RunUMAP(object.inter, reduction = "pca", dims = 1:30)

#Visualize
object.inter$border.zones <- object.inter$final_annotations_w_activity
DimPlot(object.inter , group.by = 'border.zones', cols= c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Interface' = 'grey', 
                                                                          "High" = '#FF3C31', "Medium" = 'pink',"Low" = "#5B8AFD"))

#Save object interface 
saveRDS(object.inter, 'Seurat_objects/object_inter.rds')

###### 3. Global 2-state annotation ########################################

#### 3.1 Find Optimal resolution
object.inter <- FindClusters(object.inter, verbose = FALSE, resolution = 0.05)
DimPlot(object.inter, group.by = 'integrated_snn_res.0.05' , cols = c('1' = 'red', '0' = '#28B0B0'))

# Calculate silhoutte 
embeddings <- Embeddings(object = object.inter, reduction = "pca")[, 1:30]
clusters <- object.inter$integrated_snn_res.0.05
dist_matrix <- dist(embeddings)
sil <- silhouette(as.numeric(clusters), dist_matrix)

#Plot silhoutte
plot(sil, col = 1:length(unique(clusters)), border = NA, main = "Silhouette Plot")

#### 3.2 Annotate global 2-state model 

object.inter$global_two_state <- 'Active'
object.inter$global_two_state[object.inter$integrated_snn_res.0.05 == '0'] <- 'Passive'
DimPlot(object.inter, group.by = 'global_two_state', cols = c('Active' = 'red', 'Passive' = '#3C50B1'))

#Plot proportion of spots in each cluster
df.to.plot <- object.inter@meta.data[,c('orig.ident', 'global_two_state')]
colnames(df.to.plot)[1] <- 'Sample'
df.to.plot$Sample <-gsub('sample_', '', df.to.plot$Sample)

#stacked barplot
ggplot(df.to.plot, aes(x = Sample, fill = global_two_state)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('Active' = 'red', 'Passive' = '#3C50B1')) +
  labs(
    x = "Sample",
    y = "Proportion",
    fill = "Activity State",
    title = "Proportion of States by Sample"
  ) +
  theme_minimal()

#normal barplot
ggplot(df.to.plot, aes(x = Sample, fill = global_two_state)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c('Active' = 'red', 'Passive' = '#3C50B1')) +
  labs(
    x = "Sample",
    y = "Spots",
    fill = "Activity State",
    title = "TSLI Spots per States by Sample"
  ) +
  theme_minimal()


#### 3.3 Differential Expression and Elevated Pathways by Cluster

#Find markers
DefaultAssay(object.inter) <- "Spatial"
object.inter <- NormalizeData(object.inter)
Idents(object.inter) <- "global_two_state"

# Run DE using with sample batch adjustment
de_markers <- FindMarkers(
  object = object.inter,
  ident.1 = "Active",
  ident.2 = "Passive",
  latent.vars = "orig.ident", 
  test.use = "LR",          
  min.pct = 0.1,
  logfc.threshold = 0.25
)

#Get volcano plot
res_tableOE_tb <- de_markers %>% 
  mutate(threshold_OE =p_val_adj < 0.05 & abs(avg_log2FC) >= 0.58) %>%
  mutate(
    p_val_adj_clean = ifelse(p_val_adj == 0, .Machine$double.xmin, p_val_adj),
    neg_log10_padj  = -log10(p_val_adj_clean)
  )
ggplot(res_tableOE_tb) +
  geom_point(aes(x = avg_log2FC, y = neg_log10_padj, colour = threshold_OE)) +
  ggtitle("DE Genes Active vs. Passive Activity") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))



#GO analysis of DEs
cluster.markers <- res_tableOE_tb %>%filter(p_val_adj < 0.05 & avg_log2FC > 0.58)
ego <- enrichGO(gene = rownames(cluster.markers) , 
                universe = Features(object.inter),
                keyType = "SYMBOL",
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
cluster_summary <- data.frame(ego)
## Dotplot 
dotplot(ego, showCategory=15, font.size = 12, title = 'TOP 15 GO OE in Passive Cluster ')
write.csv(as.data.frame(cluster_summary), '/path_to_data/ActiveStateDEs_GO.csv', 
          quote = F)

#Get upregulated Hallmakrs
gene_sets_df <- msigdbr(species = 'Homo sapiens', collection = 'H')
gene_sets <- gene_sets_df %>%
  dplyr::select(gs_name, gene_symbol)
x <- enricher(rownames(cluster.markers), TERM2GENE = gene_sets, universe = Features(object.inter))
cluster_summary <- data.frame(x)
barplot(x , showCategory = 15, title = paste("Upregulated Cancer Hallmarks (Passive Cluster)"))
write.csv(as.data.frame(cluster_summary), '/path_to_data/ActiveStateDEs_enriched_pathways.csv', 
          quote = F)

#save DE results
write.csv(as.data.frame(res_tableOE_tb), '/path_to_data/ActiveStateDEs.csv', 
          quote = F)
res_tableOE_tb <- read.csv('/path_to_data/ActiveStateDEs.csv')
rownames(res_tableOE_tb) <- res_tableOE_tb$X



###### 4. Find tri-state classification #####################################
object.inter <- FindClusters(object.inter, verbose = FALSE, resolution = 0.2)
DimPlot(object.inter , group.by = "integrated_snn_res.0.2")

#Make barplot of 'connecting bridge' 
bridge.spots <- as.data.frame(table(object.inter$integrated_snn_res.0.8, object.inter$orig.ident)[8,])
colnames(bridge.spots) <- 'Spots'
rownames(bridge.spots) <- letters
bridge.spots %>%
  rownames_to_column(var = "Sample") %>%
  ggplot(aes(x = reorder(Sample, -Spots), y = Spots, fill = Sample)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Sample", y = "Spots", title = 'Spots in Conecting Bridge')

#conecting bridge with resolution 0.2
object.inter$global_tri_state <- 'Active'
object.inter$global_tri_state[object.inter$integrated_snn_res.0.2 == 0] <- 'Passive'
object.inter$global_tri_state[object.inter$integrated_snn_res.0.2 %in% c(2,3,4)] <- 'Intermediate'
DimPlot(object.inter , group.by = "global_tri_state")


object.inter$global_tri_state_v2 <- 'Active'
object.inter$global_tri_state_v2[object.inter$integrated_snn_res.0.2 %in% c(0,2,4)] <- 'Passive'
object.inter$global_tri_state_v2[object.inter$integrated_snn_res.0.2 %in% c(3)] <- 'Intermediate'
DimPlot(object.inter , group.by = "global_tri_state_v2")

saveRDS(object.inter, '/path_to_data/global_3_state/object_inter.rds')

#Make barplot of bridge spots
bridge.spots <- as.data.frame(table(object.inter$integrated_snn_res.0.2, object.inter$orig.ident)[3,])
colnames(bridge.spots) <- 'Spots'
rownames(bridge.spots) <- letters
bridge.spots %>%
  rownames_to_column(var = "Sample") %>%
  ggplot(aes(x = reorder(Sample, -Spots), y = Spots, fill = Sample)) +
  geom_col() +
  theme_minimal() +
  labs(x = "Sample", y = "Spots", title = 'Spots in Conecting Bridge')


#Subset and visulize by sample
object.inter.subset <- subset(object.inter, subset = orig.ident == 'sample_E')
DimPlot(object.inter.subset , group.by = "border.zones")
DimPlot(object.inter.subset , group.by = "integrated_snn_res.0.8")
DimPlot(object.inter.subset , group.by = "global_tri_state")

table(object.inter.subset$global_tri_state, object.inter.subset$border.zones)['Intermediate',]
table(object.inter.subset$integrated_snn_res.0.8, object.inter.subset$border.zones)['7',]
SpatialPlot(object.inter.subset, group.by = 'global_tri_state', pt.size.factor = 5500, image.alpha = 0.5,
            cols = cols)
SpatialPlot(object.inter.subset, group.by = 'border.zones', pt.size.factor = 5500, image.alpha = 0.5,
            cols = cols)

res.cols <- rep('red', length(unique(object.inter.subset$integrated_snn_res.0.8)))
res.cols[8]


###### 5. Calculate Activity Score #####################################

# Calculate Activity Scores using only 'Active' DEs
active.DEs <- rownames(cluster.markers)
DefaultAssay(object.inter) <- "Spatial" 
object.inter <- NormalizeData(object.inter)

#Run AddModuleScore
object.inter <- AddModuleScore(
  object = object.inter,
  features = list(active.score = active.DEs),
  name = "Active_Score", 
  assay = "Spatial",
  search = TRUE          
)

#Visualize
VlnPlot(object.inter, 'Active_Score1', group.by = 'global_two_state', pt.size = 0, cols = c('Active' = 'red', 'Passive' = '#3C50B1'))
FeaturePlot(object.inter, 'Active_Score1')

#Save activity annotation
write.csv(as.data.frame(object.inter$Active_Score1), '/path_to_data/activity_annotation.csv', 
          quote = F)


###### 6. Make plots for every sample #####################################

#Spatial plot of every sample 
for(sample.num in 1:9){
  #subset specific sample
  sample <- letters[sample.num]
  object.inter.subset <- subset(object.inter, subset = orig.ident == paste0('sample_', sample))
  
  #save global two-state plots
  out.dir <- paste0('/path_to_data/global_2_state/', sample)
  dir.create(out.dir, showWarnings = F)
  SpatialPlot(object.inter.subset, group.by = 'global_two_state', pt.size.factor = size[sample.num],
              image.alpha = 0.5, cols = cols)
  ggsave(paste0(out.dir, '/two_state_spatial_plot.png'), units = 'in', width = 5, height = 7)
  VlnPlot(object.inter.subset, features = 'activity.score', group.by = 'global_two_state',
          cols = cols)
  ggsave(paste0(out.dir, '/activity_score_violin_plot.png'), units = 'in', width = 7, height = 7)
  
  #save global three-state plots
  out.dir <- paste0('/path_to_data/global_3_state/', sample)
  dir.create(out.dir, showWarnings = F)
  SpatialPlot(object.inter.subset, group.by = 'global_tri_state', pt.size.factor = size[sample.num],
              image.alpha = 0.5, cols = cols)
  ggsave(paste0(out.dir, '/tri_state_spatial_plot.png'), units = 'in', width = 5, height = 7)
  VlnPlot(object.inter.subset, features = 'activity.score', group.by = 'global_tri_state',
          cols = cols)
  ggsave(paste0(out.dir, '/activity_score_violin_plot.png'), units = 'in', width = 7, height = 7)
  
}




