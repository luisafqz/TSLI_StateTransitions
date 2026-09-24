######## Combine outhouse and inhouse into one object

library(Seurat)
library(org.Hs.eg.db)
library(DOSE)
library(annotables)
library(msigdbr)
library(clusterProfiler)


####### 1. Integrate external and internal BC dataset ##########################
#set samples
samples <- c("A", "B", "C", "D", "E", "G", "I", "K", "M", "T1", "T2", "T3")

#read all samples 
object.inter.list <- list()

for(sample in samples){
  #Load samples
  if(sample %in% c('T1', 'T2', 'T3')){
    base.dir <- '/path_to_data/GSE285715/'
    object.inter <- readRDS(paste0(base.dir, sample, '/boder_activity.rds'))
    object.inter$orig.ident <- sample 
    object.inter.list[[sample]] <- object.inter
  }else{
    object <- readRDS(paste0('Seurat_objects/object_inter.rds'))
    object.inter <- subset(object, subset = orig.ident == paste0('sample_', sample))
    object.inter.list[[sample]] <- object.inter
  }
}

# Merge objects
seurat_list.copy <- object.inter.list
seurat_list.copy[[1]] <- NULL
obj <- object.inter.list[[1]]
rm(object.inter.list)
merged_obj <- merge(x = obj , y = c(unlist(seurat_list.copy)),
                    project = "Integrated_Project")
DefaultAssay(merged_obj) <- "Spatial"


# Run SCTransform 
rm(seurat_list.copy, obj)
merged_obj <- SCTransform(merged_obj, assay = "Spatial", verbose = FALSE)

# Run PCA 
merged_obj <- RunPCA(merged_obj, assay = "SCT", verbose = FALSE)

# Integrate 
merged_obj <- IntegrateLayers(
  object = merged_obj, 
  method = HarmonyIntegration, 
  normalization.method = "SCT", 
  assay = "SCT",                 
  orig.reduction = "pca", 
  new.reduction = "integrated.dr",
  verbose = TRUE,
  k.weight = 50
)

# Find neighbors and clusters 
merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.dr", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 1)

# Run UMAP 
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.dr", dims = 1:30)

#Visualize
DimPlot(merged_obj, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
DimPlot(merged_obj, reduction = "umap", group.by = "orig.ident")


#test different clusters 
merged_obj <- FindClusters(merged_obj, reduction = "integrated.dr", dims = 1:30, resolution = )
DimPlot(merged_obj, reduction = "umap", group.by = "SCT_snn_res.0.04", label = TRUE)
DimPlot(merged_obj, reduction = "umap", group.by = "SCT_snn_res.1.5", label = TRUE)


#### 1.1 Add two state annotation and identify which cluster is 'Active' and 'Passive'

object_w_anno <- readRDS('/Users/lfq7/Documents/Lab/Projects/Border_transition/global_states/interface_object.rds')
object_w_anno_subset <- subset(object_w_anno, subset = global_two_state %in% c('Active', 'Passive'))

activity.by.cluster <- table(merged_obj$SCT_snn_res.1, merged_obj$interface_annot)[,c('Active', 'Passive')]
active.to.passive.ratio <- activity.by.cluster[,c('Active')]/(activity.by.cluster[,c('Passive')] + 1)
active.clusters <- as.numeric(names(active.to.passive.ratio[active.to.passive.ratio > 1]))
merged_obj$global_annotation <- 'Passive'
merged_obj$global_annotation[merged_obj$SCT_snn_res.1 %in% active.clusters] <- 'Active'
DimPlot(merged_obj, reduction = "umap", group.by = "global_annotation", cols = cols)

#Plot proportion of spots in each cluster
df.to.plot <- merged_obj@meta.data[,c('orig.ident', 'global_annotation')]
colnames(df.to.plot)[1] <- 'Sample'
df.to.plot$Sample <-gsub('sample_', '', df.to.plot$Sample)

#stacked barplot
ggplot(df.to.plot, aes(x = Sample, fill = global_annotation)) +
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


#make boxplot 
df.to.plot$Cohort <- 'In-house'
df.to.plot$Cohort[df.to.plot$Sample %in% c('T1', 'T2', 'T3')] <- 'Out-house'

#reshape df to get percentage of each sample
df_pct <- df.to.plot %>%
  count(Sample, Cohort, global_annotation) %>%
  tidyr::pivot_wider(
    names_from = global_annotation, 
    values_from = n, 
    values_fill = 0
  ) %>%
  mutate(pct_active = (Active / (Active + Passive)) * 100)

#plot
ggplot(df_pct, aes(x = Cohort, y = pct_active, fill = Cohort)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.8) + 
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 100)) +
  labs(
    title = "Percentage of Active TSLIs by Cohort",
    x = "Cohort",
    y = "Active TSLIs (%)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#Visualize 
merged_obj.subset <- subset(merged_obj, subset = orig.ident == 'T3')
SpatialFeaturePlot(merged_obj.subset, 'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION1', pt.size = 3, 
                   image.alpha = 0.5)


#### 1.2 Perform DE analysis of two clusters 
#Find markers
DefaultAssay(merged_obj) <- "Spatial"
merged_obj <- NormalizeData(merged_obj)
Idents(merged_obj) <- "global_annotation"

# Run DE using with sample batch adjustment
merged_obj <- JoinLayers(merged_obj)
de_markers <- FindMarkers(
  object = merged_obj,
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
dotplot(ego, showCategory=15, font.size = 12, title = 'TOP 15 GO OE in Active Cluster ')


#Get upregulated Hallmakrs
gene_sets_df <- msigdbr(species = 'Homo sapiens', collection = 'H')
gene_sets <- gene_sets_df %>%
  dplyr::select(gs_name, gene_symbol)
x <- enricher(rownames(cluster.markers), TERM2GENE = gene_sets, universe = Features(object.inter))
cluster_summary <- data.frame(x)
barplot(x , showCategory = 15, title = paste("Upregulated Cancer Hallmarks (Active Cluster)"))

#Spatial plots
object.subset <- subset(merged_obj, subset = orig.ident == 'T3')
SpatialPlot(object.subset, group.by = 'global_annotation', cols = cols, pt.size.factor = 3,
            image.alpha = 0.5)


####### 2.3 Calculate activity Score
# Calculate Activity Scores using only upregulated pathways
active.markers <- read.csv('/Users/lfq7/Documents/Lab/Projects/Border_transition/global_states/ActiveStateDEs.csv')
rownames(active.markers) <- active.markers$X 
cluster.markers <- active.markers %>%filter(p_val_adj < 0.05 & avg_log2FC > 0.58)
active.DEs <- rownames(cluster.markers)
DefaultAssay(merged_obj) <- "Spatial" 
merged_obj <- NormalizeData(merged_obj)

# 3. Run AddModuleScore
merged_obj <- AddModuleScore(
  object = merged_obj,
  features = list(active.score = active.DEs),
  name = "Active_Score", 
  assay = "Spatial",
  search = TRUE          
)

FeaturePlot(merged_obj, features = c('Active_Score1'))
VlnPlot(merged_obj, features = c('Active_Score1'), group.by = 'global_annotation', pt.size = 0, cols = cols)

#save results
saveRDS(merged_obj, '/Users/lfq7/Documents/Lab/Projects/Border_transition/Outhouse_samples/combines_outhouse_and_inhouse.rds')
merged_obj <- readRDS('/Users/lfq7/Library/CloudStorage/OneDrive-RutgersUniversity/Lab/Projects/Bladder_Cancer/Border_transition/Outhouse_samples/combines_outhouse_and_inhouse.rds')



###### 2. Repeat analysis using only outhouse samples ####################################

#### 2.1 Integrate External Samples
#set samples
samples <- c("T1", "T2", "T3")

#read all samples 
object.inter.list <- list()

for(sample in samples){
  #Load samples
  if(sample %in% c('T1', 'T2', 'T3')){
    base.dir <- '/path_to_data/GSE285715/'
    object.inter <- readRDS(paste0(base.dir, sample, '/boder_activity.rds'))
    object.inter$orig.ident <- sample 
    object.inter.list[[sample]] <- object.inter
  }else{
    object <- readRDS(paste0('Seurat_objects/object_inter.rds'))
    object.inter <- subset(object, subset = orig.ident == paste0('sample_', sample))
    object.inter.list[[sample]] <- object.inter
  }
}

# Merge objects
seurat_list.copy <- object.inter.list
seurat_list.copy[[1]] <- NULL
obj <- object.inter.list[[1]]
rm(object.inter.list)
merged_obj <- merge(x = obj , y = c(unlist(seurat_list.copy)),
                    project = "Integrated_Project")
DefaultAssay(merged_obj) <- "Spatial"


# Run SCTransform 
rm(seurat_list.copy, obj)
merged_obj <- SCTransform(merged_obj, assay = "Spatial", verbose = FALSE)

# Run PCA 
merged_obj <- RunPCA(merged_obj, assay = "SCT", verbose = FALSE)

# Integrate 
merged_obj <- IntegrateLayers(
  object = merged_obj, 
  method = HarmonyIntegration, 
  normalization.method = "SCT", 
  assay = "SCT",                 
  orig.reduction = "pca", 
  new.reduction = "integrated.dr",
  verbose = TRUE,
  k.weight = 50
)

# Find neighbors and clusters 
merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.dr", dims = 1:30)
merged_obj <- FindClusters(merged_obj, resolution = 1)

# Run UMAP 
merged_obj <- RunUMAP(merged_obj, reduction = "integrated.dr", dims = 1:30)

#Visualize
DimPlot(merged_obj, reduction = "umap", group.by = "border_activity", cols = cols)
DimPlot(merged_obj, reduction = "umap", group.by = "orig.ident")

DimPlot(merged_obj, reduction = "umap", group.by = "border_activity", cols = cols)

merged_obj <- FindClusters(merged_obj, resolution = 0.1)
DimPlot(merged_obj, reduction = "umap", group.by = "SCT_snn_res.0.1")

merged_obj$annotation <- 'Passive'
merged_obj$annotation[merged_obj$SCT_snn_res.0.1 == 1 ] <- 'Active'
DimPlot(merged_obj, reduction = "umap", group.by = "annotation", cols = cols)

# get propotion of spots in each state
df.to.plot <- merged_obj@meta.data[,c('orig.ident', 'annotation')]
colnames(df.to.plot)[1] <- 'Sample'
df.to.plot$Sample <-gsub('sample_', '', df.to.plot$Sample)

#stacked barplot
ggplot(df.to.plot, aes(x = Sample, fill = annotation)) +
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

##### 2.2 DE Analysis 
DefaultAssay(merged_obj) <- "Spatial"
merged_obj <- NormalizeData(merged_obj)
Idents(merged_obj) <- "annotation"

# Run DE using with sample batch adjustment
merged_obj <- JoinLayers(merged_obj)
de_markers <- FindMarkers(
  object = merged_obj,
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

res_tableOE_tb

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
#Dotplot 
dotplot(ego, showCategory=15, font.size = 12, title = 'TOP 15 GO OE in Active Cluster ')

#Get upregulated Hallmakrs
gene_sets_df <- msigdbr(species = 'Homo sapiens', collection = 'H')
gene_sets <- gene_sets_df %>%
  dplyr::select(gs_name, gene_symbol)
x <- enricher(rownames(cluster.markers), TERM2GENE = gene_sets, universe = Features(object.inter))
cluster_summary <- data.frame(x)
barplot(x , showCategory = 15, title = paste("Upregulated Cancer Hallmarks (Active Cluster)"))

#Spatial plots
object.subset <- subset(merged_obj, subset = orig.ident == 'T3')
SpatialPlot(object.subset, group.by = 'annotation', cols = cols, pt.size.factor = 3,
            image.alpha = 0.5)

#UMAP plots 
FeaturePlot(merged_obj, c('HALLMARK_APICAL_JUNCTION1'))

#save out-house only analysis 
saveRDS(merged_obj, '/path_to_data/combined_outhouse.rds')
merged_obj <- readRDS('/path_to_data/combined_outhouse.rds')


#### 3. Extract data for COMMOT ####################################################
combined.object <- readRDS('/path_to_data/combines_outhouse_and_inhouse.rds')

for(sample in c('T1', 'T2', 'T3')){
  #add border annotation
  base.dir <- '/path_to_data/GSE285715/'
  object.inter <- readRDS(paste0(base.dir, sample, '/Seurat_obj_w_anno.rds'))
  annotation <- data.frame(Tumor_regions = object.inter$final_annotations)
  
  #add region annotation
  combined.object.subset <- subset(combined.object, subset = orig.ident == sample)
  object.inter$border_annotation <- object.inter$final_annotations
  object.inter$border_annotation[object.inter$final_annotations == 'Border'] <- combined.object.subset$global_annotation
  
  annotation$Border_annotation <- object.inter$border_annotation
  
  write.csv(annotation, paste0(out.dir, sample, '_border_annotations.csv'))
}


