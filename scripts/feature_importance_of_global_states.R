### Get important pathways for clustering  


library(Seurat)
library(dplyr)
library(org.Hs.eg.db)
library(DOSE)
library(annotables)
library(msigdbr)
library(clusterProfiler)

source('/path_to_CALISTA/ModuleB.R')
source('/path_to_CALISTA/ModuleA.R')
source('/path_to_data/feature_importance_formulas.R')


#Set parameters
input.dir <- 'path_to_data'
letters <- c("A", "B", "C", "D", "E", "G", "I", "K", "M")
samples <- c("sampleA", "sampleB", "sampleC", "sampleD", "sampleE", "sampleG", "sampleI", 
             "sampleK", "sampleM")
size <- c(rep(10500, 4), 5500, 3000, 4000, 5500, 3500)
cols <- c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Interface' = 'grey', 
          "Active" = '#FF3C31', "Passive" = "#5B8AFD", 'Intermediate' = '#b5179e')

#read RDS 
object.inter <- readRDS('/path_to_data/global_3_state/object_inter.rds')
DimPlot(object.inter , group.by = "global_tri_state")

##### 1. Get Proportion of spots in each activity cluster ########################
#Plot proportion of spots in each cluster
df.to.plot <- object.inter@meta.data[,c('orig.ident', 'global_tri_state')]
colnames(df.to.plot)[1] <- 'Sample'
df.to.plot$Sample <-gsub('sample_', '', df.to.plot$Sample)

#stacked barplot
ggplot(df.to.plot, aes(x = Sample, fill = global_tri_state)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c('Active' = 'red', 'Passive' = '#3C50B1', 'Intermediate' = '#b5179e')) +
  labs(
    x = "Sample",
    y = "Proportion",
    fill = "Activity State",
    title = "Proportion of States by Sample"
  ) +
  theme_minimal()

#unstacked barplot
ggplot(df.to.plot, aes(x = Sample, fill = global_tri_state)) +
  geom_bar(position = "stack") +
  scale_fill_manual(values = c('Active' = 'red', 'Passive' = '#3C50B1', 'Intermediate' = '#b5179e')) +
  labs(
    x = "Sample",
    y = "Spots",
    fill = "Activity State",
    title = "TSLI Spots per States by Sample"
  ) +
  theme_minimal()

#plot only at intermediate state
ggplot(df.to.plot[df.to.plot$global_tri_state == 'Intermediate',], 
       aes(x = reorder(Sample, Sample, FUN = function(x) -length(x)), fill = Sample)) +
  geom_bar(position = "stack") +
  labs(
    x = "Sample",
    y = "Spots",
    fill = "Sample",
    title = "TSLI Spots in Intermediate State"
  ) +
  theme_minimal()

#Save annotation
write.csv(as.data.frame(object.inter$global_tri_state), '/path_to_data/tri_state_border_annotation.csv', quote = F)

#### 2. Get feature importance for classification ############################################

#get Hallmarks df 
cancer.hallmarks.meta <- object.inter[[]][,grepl('HALLMARK', colnames(object.inter@meta.data))]
normalized.data <- as.data.frame(scale(cancer.hallmarks.meta))
normalized.data$Cluster <- as.factor(object.inter$global_tri_state)
colnames(normalized.data) <- gsub('EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT', colnames(normalized.data))
colnames(normalized.data) <- gsub('HALLMARK_', '', colnames(normalized.data))

#get importances
MDA.list <- run.PIMP.MDA(normalized.data, get.important.pathways = T)
MDI.list <- run.PIMP.MDI(normalized.data, get.important.pathways = T)

#save lists 
saveRDS(MDA.list, '/path_to_data/MDA_list.rds')
saveRDS(MDI.list, '/path_to_data/MDI_list.rds')

# Get matching MDA and MDI
match.sig <- intersect(MDI.list[['v1']], MDA.list[['v1']])

#Get expression of significant Pathways
FeaturePlot(object.inter, features = match.sig)


##### 3. Find Marker Genes ############################################################################

#Get marker genes
Idents(object.inter) <- object.inter$global_tri_state
DefaultAssay(object.inter) <- "Spatial"
cluster.markers <- FindAllMarkers(object.inter, min.pct = 0.1,
                                  latent.vars = "orig.ident", 
                                  test.use = "LR",          
                                  logfc.threshold = 0.25)
top <- cluster.markers  %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

# Plot top genes Heatmap 
DefaultAssay(object.inter) <- "SCT"
object.inter <- ScaleData(object.inter, features = top$gene)
DoHeatmap(object.inter, features = top$gene) + NoLegend()


#Run GO analysis for specified group
group <- 'Passive'
res_table_OE_tb <- cluster.markers %>% filter(cluster == group ) %>%
  filter(p_val_adj < 0.05 & avg_log2FC >= 0.58)
ego <- enrichGO(gene = rownames(res_table_OE_tb),
                keyType = "SYMBOL",
                universe = Features(object.inter),
                OrgDb = org.Hs.eg.db, 
                ont = "ALL", 
                pAdjustMethod = "fdr", 
                pvalueCutoff = 0.05)


#make dotpot
dotplot(ego, showCategory=15, font.size = 12, title = paste0('TOP 15 GO for "', group, '" Marker Genes'))

#get Hallmark enrichment
x <- enricher(gsub('.1', '', rownames(res_table_OE_tb)), TERM2GENE = gene_sets, universe = Features(object.inter))
barplot(x , showCategory = 15, title = paste0("Upregulated Cancer Hallmarks (", group,")"))



##### 4. DE Analysis ############################################################################

#Find markers of specified clusters
DefaultAssay(object.inter) <- "Spatial"
object.inter <- NormalizeData(object.inter)
Idents(object.inter) <- "global_tri_state"

de_markers <- FindMarkers(
  object = object.inter,
  ident.1 = "Intermediate",
  ident.2 = "Active",
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
  ggtitle("DE Genes Intermediate vs. Passive Activity") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))
write.csv(as.data.frame(res_tableOE_tb), '/path_to_data/Intermediate_vs_Active_StateDEs.csv', 
          quote = F)

#GO Analysis
cluster.markers <- res_tableOE_tb %>%filter(p_val_adj < 0.05 & avg_log2FC > -0.58)
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
dotplot(ego, showCategory=10, font.size = 10, title = 'TOP 10 GO OE in Intermediate Cluster ')
write.csv(as.data.frame(cluster_summary), '/path_to_data/Intermediate_vs_Active_StateDEs_GOs.csv', 
          quote = F)

#Get upregulated Hallmakrs
gene_sets_df <- msigdbr(species = 'Homo sapiens', collection = 'H')
gene_sets <- gene_sets_df %>%
  dplyr::select(gs_name, gene_symbol)
x <- enricher(rownames(cluster.markers), TERM2GENE = gene_sets, universe = Features(object.inter))
cluster_summary <- data.frame(x)
barplot(x , showCategory = 15, title = paste("Upregulated Cancer Hallmarks (Intermediate Cluster)"))
write.csv(as.data.frame(cluster_summary), '/path_to_data/Intermediate_vs_Active_StateDEs_enriched_hallmarks.csv', 
          quote = F)



##### 5. Get Silhouette statistic ###########################################

#Get embeddings 
embeddings <- Embeddings(object = object.inter, reduction = "pca")[, 1:30]

#Get clusters
DimPlot(object.inter, group.by = 'integrated_snn_res.0.2' )
clusters <- as.numeric(as.factor(object.inter$integrated_snn_res.0.2))
dist_matrix <- dist(embeddings)
sil <- silhouette(as.numeric(clusters), dist_matrix)

#Plot silhoutte
plot(sil, col = 1:length(unique(clusters)), border = NA, main = "Resolution 0.2")



#### 6. Get feature importance for classification of two-state ############################################

#extract object
object <- readRDS('/path_to_data/interface_object.rds')
object.inter <- subset(object, subset = interface_w_border == 'Border')

#get df 
cancer.hallmarks.meta <- object.inter[[]][,grepl('HALLMARK', colnames(object.inter@meta.data))]
normalized.data <- as.data.frame(scale(cancer.hallmarks.meta))
normalized.data$Cluster <- as.factor(object.inter$global_two_state)
colnames(normalized.data) <- gsub('EPITHELIAL_MESENCHYMAL_TRANSITION', 'EMT', colnames(normalized.data))
colnames(normalized.data) <- gsub('HALLMARK_', '', colnames(normalized.data))

#get importances
MDA.list <- run.PIMP.MDA(normalized.data, get.important.pathways = T)
MDI.list <- run.PIMP.MDI(normalized.data, get.important.pathways = T)

# Get matching MDA and MDI
match.sig <- intersect(MDI.list, MDA.list)
match.sig <- gsub('EMT', 'EPITHELIAL_MESENCHYMAL_TRANSITION', match.sig)
match.sig <- paste0('HALLMARK_', match.sig)

#Get expression of significant Pathways
VlnPlot(object.inter, features = match.sig, group.by = 'global_two_state', pt.size = 0)
FeaturePlot(object.inter, features = match.sig)


