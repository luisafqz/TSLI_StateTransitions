### Statistics about local borders 

library(Seurat)
library(sp)
library(spdep)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

#set params 
letters <- c("A", "B", "C", "D", "E", "G", "I", "K", "M")
samples <- c("sampleA", "sampleB", "sampleC", "sampleD", "sampleE", "sampleG", "sampleI", 
             "sampleK", "sampleM")
size <- c(rep(10500, 4), 5500, 3000, 4000, 5500, 3500)

##### 1. Get Local Interfaces ###########################################

## 1.1 Create Neighborhood object 
#load interface object
#read object
object.interface <- readRDS('/Users/lfq7/Documents/Lab/Projects/Border_transition/global_states/interface_object.rds')

#subset to keep only TSLI region 
object.inter <- subset(object.interface, subset = interface_w_border == 'Border')

#Create coord list 
coords.list <- list()
for(i in 1:9){
  object.inter.subset <- subset(object.inter, subset = orig.ident == paste0('sample_', letters[i]))
  coords <- GetTissueCoordinates(object.inter.subset)
  coords.list[[i]] <- coords
}

#Get Neighborhood Data for all samples
subgraphs.list <- list()
subraph.sizes.list <- list()
subgraph.IDs.list <- c()
for(i in 1:9){
  #create neighborhood object
  sp.tme <- coords.list[[i]]
  colnames(sp.tme) <- c('row', 'col')
  pt.dist <- min(dist(coords.list[[i]]))
  sp::coordinates(sp.tme) <- ~ row  + col
  neib <- spdep::dnearneigh(as.matrix(sp::coordinates(sp.tme)), d1 = 1, d2 = pt.dist + pt.dist/2)
  
  #extract data 
  comp.info <- n.comp.nb(neib)
  subgraphs <- comp.info$nc
  subgraph.sizes <- table(comp.info$comp.id)
  subgraph.IDs <-as.factor(comp.info$comp.id)
  
  #Add to list 
  subgraphs.list[[i]] <- subgraphs
  subraph.sizes.list[[i]] <- subgraph.sizes
  subgraph.IDs.list <- c(subgraph.IDs.list, subgraph.IDs)
}

#Add subgraph list to object
object.inter$subgraph_ID <- subgraph.IDs.list

#Read color list 
vibrant_50 <- read.csv('/path_to_data/TSLI_colors.csv')
vibrant_50 <- vibrant_50[,2]
names(vibrant_50) <- 1:47

#Visualize and save results for all samples
sample.num <- 8
for(sample.num in 1:9){
  sample <- letters[sample.num]
  object.inter$TSLI <- object.inter$subgraph_ID
  SpatialPlot(object.inter, group.by = 'TSLI', images = samples[sample.num], 
              pt.size.factor = size[sample.num], image.alpha = 0.3, cols = vibrant_50)
  ggsave(paste0('/path_to_data/', sample, '/TSLI_annotation.png'))
}

#Visualize a specific TSLI
specific.tsli <- rep('darkgrey', 47)
names(specific.tsli) <- 1:47
specific.tsli[c(11, 8, 24)] <- vibrant_50[c(11, 8, 24)]
sample.num <- 8
sample <- letters[sample.num]
object.inter$TSLI <- object.inter$subgraph_ID
SpatialPlot(object.inter, group.by = 'TSLI', images = samples[sample.num], 
            pt.size.factor = size[sample.num], image.alpha = 0.3, cols = specific.tsli)



#### 1.2  Get subgraph statistics  
#Get subgraph plots 
subgraphs.list <- unlist(subgraphs.list)
names(subgraphs.list) <- letters
barplot(subgraphs.list, main = "TSLI's per sample", ylab = 'Count', 
        col = c('#F9766D', '#C19746', '#9AA848', '#6EB55A', '#72BCA2', '#6DB3DA', '#789AF1', '#C17CED', '#ED6CC0'))


### 1.3 Get cell variance across all samples 
cell_var <- list()

for(i in 1:9){
  #extract data
  object.inter.subset <- subset(object.inter, subset = orig.ident == paste0('sample_', letters[i]))
  h.df <- object.inter.subset@meta.data[,grep('HALLMARK', colnames(object.inter.subset@meta.data))]
  h.df$border.zones <- object.inter.subset$global_tri_state
  h.df$border.zones <- factor(h.df$border.zones)
  
  #Create hallmarks data frame 
  var.df <- as_tibble(h.df) %>% group_by(border.zones) %>% summarise(across(where(is.numeric), var, na.rm = TRUE))
  
  #Get border zone variance
  h.df$cell_variance <- apply(h.df %>% select(where(is.numeric)), 1, var)
  
  # Then find the average cell variance per ID
  grouped_cell_var <- h.df %>%
    group_by(border.zones) %>%
    summarise(mean_cell_variance = mean(cell_variance))
  
  cell_var[[i]] <- grouped_cell_var
}

cell.var.mat <- c()
for(i in 1:9){
  cell.var.mat <- cbind(cell.var.mat,  as.numeric(cell_var[[i]]$mean_cell_variance))
}
rownames(cell.var.mat) <- c('Active', 'Intermediate', 'Passive')
colnames(cell.var.mat) <- letters


#Plot Heatmap
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
Heatmap(cell.var.mat, name = 'Mean var', col = col_fun, cluster_columns = F, 
        cluster_rows = F, column_title = 'Mean Cell Variance of Activity State')
write.csv(cell.var.mat, '/path_to_data/mean_cell_variance_of_activity_state.csv')


##### 2. Correlation between border spots and tissue size  #######################################

#extract data frame
object <- readRDS('Seurat_objects/object_with_annotations.rds')
df.to.plot <- object@meta.data[,c("orig.ident", 'final_annotations')]
to.plot <- as.data.frame(table(df.to.plot$orig.ident, df.to.plot$final_annotations))
colnames(to.plot) <- c('Sample', 'Region', 'Count')
to.plot[to.plot$Var2 != "Border",]

#plot Composition of Regions
ggplot(to.plot, aes(x = Sample, y = Count, fill = Region)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Border' = 'black')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt labels for readability
  labs(title = "Region Composition per Sample",
       x = "Sample ID",
       y = "Number of Spots")

ggplot(to.plot, aes(x = Sample, y = Count, fill = Region)) +
  geom_bar(stat = "identity", position = 'fill') +
  scale_fill_manual(values = c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Border' = 'black')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt labels for readability
  labs(title = "Region Composition per Sample",
       x = "Sample ID",
       y = "Number of Spots")

#Plot tumor to stroma ratio
ggplot(to.plot[to.plot$Region != "Border",], aes(x = Sample, y = Count, fill = Region)) +
  geom_bar(stat = "identity", position = 'fill') +
  scale_fill_manual(values = c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Border' = 'black')) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + # Tilt labels for readability
  labs(title = "Region Composition per Sample",
       x = "Sample ID",
       y = "Tumor-Stroma Ratio")

#Get correlation of border spots per sample and border spots
to.plot <- cbind(table(df.to.plot$orig.ident),table(df.to.plot$orig.ident, df.to.plot$final_annotations)[,1]) 
colnames(to.plot) <- c("Spots.in.Sample", "Border.Spots")
to.plot <- as.data.frame(to.plot)

#plot spots vs border spots 
plot(to.plot[,1], to.plot[,2], pch = 19, cex = 1, ylab = 'Border Spots', xlab = 'Spots in Sample')

#fit and add trend line
fit <- lm(Border.Spots ~ Spots.in.Sample, data = to.plot)
abline(fit, col = "red", lwd = 2)

#Calculate correlation and add to plot
r_value <- cor(to.plot[,1], to.plot[,2], use = "complete.obs")
text(x = min(to.plot[,1]), y = max(to.plot[,2]) -5, 
     labels = paste("R =", round(r_value, 2)), 
     pos = 4, cex = 1.2, col = "red")



##### 3. Homogeneity across all samples  #######################################

#get variance by region 
h.df <- object.inter@meta.data[,grep('HALLMARK', colnames(object.inter@meta.data))]
h.df$orig.ident <- object.inter$orig.ident
h.df$orig.ident <- factor(h.df$orig.ident)

#Get border zone variance
h.df$cell_variance <- apply(h.df %>% select(where(is.numeric)), 1, var)

# Then find the average cell variance per ID
grouped_cell_var <- h.df %>%
  group_by(orig.ident) %>%
  summarise(mean_cell_variance = mean(cell_variance))

#Get barplot 
variance.list <- grouped_cell_var$mean_cell_variance
names(variance.list ) <- letters
barplot(variance.list, main = "Border Spots Mean Variance per Sample", ylab = 'Variance', 
        col = c('#F9766D', '#C19746', '#9AA848', '#6EB55A', '#72BCA2', '#6DB3DA', '#789AF1', '#C17CED', '#ED6CC0'))


##### 4. Global activity Score  #######################################

source('/path_to_data/ModuleA.R')

#add activity score
object.inter$activity.score <- object.inter$Active_Score1

#plot active score for all samples
for(i in 1:9){
  sample <- letters[i]
  out.dir <- paste0(out.dir <- '/Users/lfq7/Documents/Lab/Projects/Border_transition/global_states/global_2_state/', sample, '/active_score_spatial_plot.png')
  object.inter.subset <- subset(object.inter, subset = orig.ident == paste0('sample_', letters[i])) 
  SpatialFeaturePlot(object.inter.subset, 'Active_Score1', pt.size.factor = size[i]) 
  ggsave(out.dir)
}

### 4.1 Calculate autocorrelation of activity scores 
activity.scores.list <- list()
autocorrelation.list <- list()
tsli.lengths.list <- list()

#loop through every sample
for(sample.num in 1:9){
  #Subset object
  sample <- letters[sample.num]
  object <- subset(object.inter, subset = orig.ident == paste0('sample_', sample))
  
  #Get sizes of TSLI and subset to keep only those greater than 5
  sizes <- table(object$subgraph_ID)
  TSLI.to.keep <- names(sizes[sizes > 10])
  object.sub <- subset(object, subset = subgraph_ID %in% TSLI.to.keep )
  
  #save length of TSLIs
  tsli.lengths.list[[sample]] <- sizes[TSLI.to.keep]
  
  #Get TSLI data
  coords.list <- list()
  x.list <- list()
  for(i in TSLI.to.keep){
    #extract TSLI object
    object.TSLI <- subset(object.sub, subset = subgraph_ID == i)
    
    #Get coordinates and annotation
    coords.list[[as.character(i)]] <- GetTissueCoordinates(object.TSLI)
    
    #Get one-hot encoded of activity
    clusters <- as.factor(object.TSLI$global_tri_state)
    if(length(levels(clusters)) > 1){
      X.mat <- model.matrix(~ clusters  - 1)
      colnames(X.mat) <- gsub('clusters', '', colnames(X.mat))
      X.mat <- as.data.frame(X.mat) 
    }else{
      X.mat <- data.frame(rep(1,length(clusters)))
      colnames(X.mat) <- paste0('clusters', levels(clusters)[1])
    }
    
    X.mat$activity.score <- object.TSLI$activity.score
    x.list[[as.character(i)]] <- X.mat
    
  }
  
  #Calculate moran
  moran.res <- list()
  for(i in TSLI.to.keep){
    if(dim(coords.list[[as.character(i)]])[1] > 10){
      scores.df <- x.list[[as.character(i)]]
      #Create a neighborhood object and extract neighbors list 
      sp.tme <- coords.list[[as.character(i)]]
      colnames(sp.tme) <- c('row', 'col')
      pt.dist <- min(dist(coords.list[[as.character(i)]]))
      sp::coordinates(sp.tme) <- ~ row  + col
      neib <- spdep::dnearneigh(as.matrix(sp::coordinates(sp.tme)), d1 = 1, d2 = pt.dist + pt.dist/2)
      neib.listw <- spdep::nb2listw(neib,  zero.policy = TRUE)
      
      #create matrix to store results
      mat <- matrix(nrow = length(colnames(scores.df)), ncol = 2)
      colnames(mat) <- c("Moran I", "p-value")
      rownames(mat) <- colnames(scores.df)
      
      #Calculate autocorrelation
      for (j in 1:length(colnames(scores.df))){
        cell.tme=scores.df[,j]
        mt <- adespatial::moran.randtest(cell.tme, neib.listw, nrepet = 9999)
        mat[j, 1] <- mt$obs[[1]]
        mat[j, 2] <- mt$pvalue
      }
      moran.res[[as.character(i)]] <- mat
    }
  }
  
  #save moran results
  autocorrelation.list[[sample]] <- moran.res
  
  #Get mean Score for each TSLI
  mean.df <- c()
  for(i in TSLI.to.keep){
    mean.res <- data.frame(TSLI.num = i, activity.score = mean(x.list[[i]]$activity.score))
    mean.df <- rbind(mean.df, mean.res)
  }
  
  #Add to list
  activity.scores.list[[sample]] <- mean.df
  
}

#format results
plot.df <- c()
for(sample.num in 1:9){
  sample <- letters[sample.num]
  plot.df <- rbind(plot.df,cbind(rep(sample, dim(activity.scores.list[[sample]])[1]), 
                                 activity.scores.list[[sample]], as.numeric(tsli.lengths.list[[sample]])))
}
colnames(plot.df)[1] <- 'Sample'
colnames(plot.df)[4] <- 'Length'


#Make scatter plot
ggplot(plot.df, aes(x = Length, y = activity.score, color = Sample)) +
  geom_point(alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    x = "Length",
    y = "Active Score",
    color = "Sample",
    title = "Activity Score vs. Length by Sample"
  )


#Format autocor df 
autocor.df <- c()
for(sample.num in 1:9){
  sample <- letters[sample.num]
  for(i in names(tsli.lengths.list[[sample]])){
    mat <- as.data.frame(autocorrelation.list[[sample.num]][[i]])
    mat <- mat[c('High', 'activity.score'),]
    mat$Score <- c('High', 'activity.score')
    mat$Sample <- sample
    mat$TSLI <- i
    mat$Length <- tsli.lengths.list[[sample]][i]
    autocor.df <- rbind(autocor.df, mat)
  }
}

#plot autocor
autocor.df %>%
  filter(Score == "activity.score") %>%
  ggplot(aes(x = Length, y = `Moran I`, color = Sample)) +
  geom_point(alpha = 0.7, size = 4) +
  theme_minimal() +
  labs(
    x = "Length",
    y = "Moran's I",
    color = "Sample",
    title = "Moran's I vs. Length of TSLI"
  )


#Add score to autocor df and compare res
score.df <- autocor.df %>%
  filter(Score == "activity.score")
score.df$Score <- plot.df$activity.score
score.df %>% ggplot(aes(x = Score, y = `Moran I`, color = Sample)) +
  geom_point(alpha = 0.7, size = 4) +
  theme_minimal() +
  labs(
    x = "Active Score",
    y = "Moran's I",
    color = "Sample",
    title = "Moran's I vs. Active Score of TSLI"
  )

#plot length vs score
score.df %>% ggplot(aes(x = Length, y = Score, color = Sample)) +
  geom_point(alpha = 0.7, size = 4) +
  theme_minimal() +
  labs(
    x = "Length",
    y = "Active Score",
    color = "Sample",
    title = "Length vs. Active Score of TSLI"
  )

#Visialize results for each sample
score.df  %>%
  filter(Sample == "B") %>%
  ggplot(aes(x = Score, y = `Moran I`, color = TSLI)) +
  geom_point(alpha = 0.7, size = 5) +
  theme_minimal() +
  labs(
    x = "Score",
    y = "Moran's I",
    color = "TSLI",
    title = "Moran's I vs. Activity Score of TSLI"
  ) + scale_color_manual(values = vibrant_50)


# Get correlation for TSLI features
cor_matrix <- score.df %>%
  group_by(Sample) %>%
  summarise(
    # Compute pairwise correlations
    cor_Moran_Score  = cor(`Moran I`, Score, use = "pairwise.complete.obs"),
    cor_Moran_Length = cor(`Moran I`, Length, use = "pairwise.complete.obs"),
    cor_Score_Length = cor(Score, Length, use = "pairwise.complete.obs")
  ) %>%
  column_to_rownames(var = "Sample") %>%
  as.matrix()

# Plot heatmap
colnames(cor_matrix) <- gsub('cor_', '',colnames(cor_matrix))
Heatmap(t(cor_matrix), name = 'Cor')

#Calculate Aggressivenes (ie. Active Margin Score)
score.df$Aggressiveness <- score.df$`Moran I` * score.df$Score

#Visualize by sample
sample <- 'C'
score.df %>%
  filter(Sample == sample) %>%
  mutate(TSLI = reorder(as.character(TSLI), Aggressiveness)) %>%
  ggplot(aes(x = TSLI, y = Aggressiveness, fill = Aggressiveness)) +
  geom_col(show.legend = FALSE) +
  scale_fill_gradient(low = "#5B8AFD", high = "#FF3C31") +
  theme_minimal(base_size = 14) +
  labs(
    x = "TSLI",
    y = "Active Margin Score",
    title = paste0("Sample ", sample, ": TSLIs Sorted by Active Margin Score")
  ) +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )




### 4.2 Make plot with combined LR and TSLI Margin Score
sample_target <- "M" 

#Extract sample data
sample_scores <- score.df %>%
  filter(Sample == sample_target) %>%
  mutate(TSLI = reorder(as.character(TSLI), Aggressiveness))

#Get TSLI Order
tsli_order <- levels(sample_scores$TSLI)

#Make Barplot
p_bar <- ggplot(sample_scores, aes(x = TSLI, y = Aggressiveness, fill = Aggressiveness)) +
  geom_col(show.legend = FALSE) +
  scale_fill_gradient(low = "#5B8AFD", high = "#FF3C31") +
  theme_minimal(base_size = 12) +
  labs(
    x = NULL, 
    y = "Active Margin Score",
    title = paste0("Sample ", sample_target, ": TSLIs Sorted by Active Margin Score")
  ) +
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 20),
    axis.title.y = element_text(size = 17)
  )

#Get TSLI colors 
built_p <- ggplot_build(p_bar)
bar_colors <- built_p$data[[1]] %>%
  select(x, fill, y) %>%
  mutate(TSLI = sample_scores$TSLI) %>% 
  select(c('TSLI', 'fill'))
write.csv(bar_colors, paste0('/path_to_data/', sample_target, 
                 '/TSLI_color_list.csv'))

#read LRs upregulated in TSLIs
heatmap.df <- read.csv(paste0('/path_to_data/', sample_target, 
                                '/upregulated_LR_per_TSLI.csv'))
heatmap.df <- heatmap.df[,-1]
colnames(heatmap.df) <- c('TSLI', 'LR')
heatmap_df <- heatmap.df %>%
  group_by(TSLI) %>% 
  pivot_longer(cols = -TSLI, names_to = "Feature", values_to = "Upreg. LR") %>%
  mutate(TSLI = factor(as.character(TSLI), levels = tsli_order)) 

#make heatmap
p_heat <- ggplot(heatmap_df, aes(x = TSLI, y = Feature, fill = `Upreg. LR`)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient(low = "#D3D3D3", high = "black", limits = c(0, NA)) +
  theme_minimal(base_size = 20) +
  labs(x = "TSLI", y = NULL, fill = "Number of Upregulated LR") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

#Format and combine both plots
p_bar <- p_bar + theme(
  plot.margin = margin(b = -10, unit = "pt")
)
p_heat <- p_heat + theme(
  plot.margin = margin(t = -10, unit = "pt")
)
p_bar / p_heat + 
  plot_layout(heights = c(8, 1)) & 
  theme(plot.margin = margin(0, 0, 0, 0)) &
  plot_annotation() & 
  theme(panel.spacing = unit(0, "lines"))
ggsave(paste0('/Users/lfq7/Documents/Lab/Projects/Border_transition/global_states/global_2_state/', sample_target, 
              '/TSLI_Active_Margin_Score_barplot.png'), units = 'in', height = 6, width = 10)


#Get sd of Aggressiveness
score.df %>%
  group_by(Sample) %>%
  summarise(
    mean_aggressiveness = mean(Aggressiveness, na.rm = TRUE),
    sd_aggressiveness   = sd(Aggressiveness, na.rm = TRUE),
    count               = n()
  )

#Plot Aggressiveness by Sample
score.df %>%
  ggplot(aes(x = Sample, y = Aggressiveness, fill = Sample)) +
  geom_boxplot(show.legend = FALSE, alpha = 0.7, outlier.colour = "red") +
  scale_fill_manual(values = unname(vibrant_50)) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Sample",
    y = "Aggressiveness",
    title = "Distribution of Aggressiveness by Sample"
  ) +
  theme(
    panel.grid.major.x = element_blank())

#Plot sd of Aggressiveness by Sample
score.df %>%
  group_by(Sample) %>%
  summarise(sd_aggressiveness = sd(Aggressiveness, na.rm = TRUE)) %>%
  ggplot(aes(x = Sample, y = sd_aggressiveness, fill = Sample)) +
  geom_col(show.legend = FALSE) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Sample",
    y = "Active Margin Score (SD)",
    title = "Active Margin Score Standard Deviation by Sample"
  ) +
  theme(
    panel.grid.major.x = element_blank()
  )

#Save TSLI annotation for Commot 
write.csv(data.frame(sample = object.inter$orig.ident, TSLI = object.inter$subgraph_ID), 
          '/path_to_data/TSLI_annotation.csv',
          quote = F)


## 4.3 Compare Cell Cycle Scores
#Get cell cycle genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

#Score cell cycle 
object.inter  <- CellCycleScoring(object.inter , s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

cell.cycle.df <- object.inter@meta.data[,c("orig.ident","subgraph_ID","S.Score", "G2M.Score")]
cell.cycle.df$orig.ident <- gsub('sample_', '', cell.cycle.df$orig.ident)
colnames(cell.cycle.df)[1:2] <- c('Sample', 'TSLI')

#calculate mean of each TLSI
mean_cell_cycle_df <- cell.cycle.df %>%
  group_by(Sample, TSLI) %>%
  summarise(
    mean_S_Score = mean(S.Score, na.rm = TRUE),
    mean_G2M_Score = mean(G2M.Score, na.rm = TRUE),
    .groups = "drop"
  )

#Add to score.df
score.df <- score.df %>% mutate(TSLI = as.character(TSLI))
mean_cell_cycle_df <- mean_cell_cycle_df %>% mutate(TSLI = as.character(TSLI))
score.df <- score.df %>%
  left_join(
    mean_cell_cycle_df %>% select(Sample, TSLI, mean_S_Score, mean_G2M_Score),
    by = c("Sample", "TSLI")
  )


#Calculate correlation
cor_matrix <- score.df %>%
  group_by(Sample) %>%
  summarise(
    cor_Active_Score_S_phase  = cor(Score, mean_S_Score, use = "pairwise.complete.obs"),
    cor_Active_Score_G2M_phase  = cor(Score, mean_G2M_Score, use = "pairwise.complete.obs"),
    cor_Active_Margin_Score_S_phase  = cor(Aggressiveness, mean_S_Score, use = "pairwise.complete.obs"),
    cor_Active_Margin_Score_G2M_phase  = cor(Aggressiveness, mean_G2M_Score, use = "pairwise.complete.obs")
    
  ) %>%
  # Convert the 'Sample' column into matrix row names
  column_to_rownames(var = "Sample") %>%
  as.matrix()
write_csv(as.data.frame(cor_matrix), '/Users/lfq7/Documents/Lab/Projects/Border_transition/supplementary_tables/cell_cycle_correaltion_res.csv')

#Plot correlation
colnames(cor_matrix) <- gsub('cor_', '',colnames(cor_matrix))
Heatmap(t(cor_matrix[,1:2]), name = 'Cor')
