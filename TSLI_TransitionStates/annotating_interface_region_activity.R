#Annotate Interface Region 

library(Seurat)
library(h2o)
library(sp)
library(gstat)
library(ggpubr)
library(ggplot2)
library(reshape2)

source('path_to_CALISTA/ModuleA.R')

#set parameters
base.dir <- '/path_to_data/'
letters <- c("A", "B", "C", "D", "E", "G", "I", "K", "M")
sample.num <- 9
samples <- c("sampleA", "sampleB", "sampleC", "sampleD", "sampleE", "sampleG", "sampleI", 
             "sampleK", "sampleM")
size <- c(rep(10500, 4), 5500, 3000, 4000, 5500, 3500)
cols <- c('Tumor' = '#F5DE95', 'Stroma' = "#DABCDF", 'Interface' = 'grey', 
          "Active" = '#FF3C31', "Int." = 'pink',"Passive" = "#5B8AFD", 'Intermediate' = '#b5179e',)

###### 1. Predict Interface annotation ################################################
#Read interface object
object.interface <- readRDS('/path_to_data/interface_object.rds')

#subset by sample
sample.num <- 9
sample <- letters[sample.num]
object.subset <- subset(object.interface , subset = orig.ident == paste0('sample_', sample))
out.dir <- paste0(base.dir, sample, '/interface_res/')
dir.create(out.dir)

#Predict using all samples 
object.subset <- object.interface 
out.dir <- paste0(base.dir, '/interface_res/')
dir.create(out.dir)

#Start h2o session
min_port = 50000
max_port = 60000
random_port <- sample(min_port:max_port, 1)
h2o::h2o.no_progress()  
h2o::h2o.init(port = random_port, bind_to_localhost = TRUE) 

#Extract input data
coords <- Seurat::GetTissueCoordinates(object.subset)[,1:2]
colnames(coords) <- c('imagerow', 'imagecol')
object.subset$interface_w_border[object.subset$interface_w_border != 'Interface'] <- object.subset$global_tri_state[object.subset$interface_w_border != 'Interface']

#coords 
coords <- c()
for (sample in letters){
  object.sample <- subset(object.interface , subset = orig.ident == paste0('sample_', sample))
  coords <- rbind(coords, Seurat::GetTissueCoordinates(object.sample )[,1:2])
}
colnames(coords) <- c('imagerow', 'imagecol')

#Set input data
input.data <- as.data.frame(cbind(Barcodes = rownames(coords), coords,
                                  object.subset@meta.data[,c('Active_Score1',colnames(object.subset@meta.data)[grepl('HALLMARK', colnames(object.subset@meta.data))])],
                                  interface_w_border = object.subset@meta.data[,c('interface_w_border')]))

#Extract data that will be use to train and test the model 
input.data <- dplyr::as_tibble(input.data)
nn.model.data <- input.data |> dplyr::filter(interface_w_border !=  'Interface')

#create h2o object
nn.model.data$interface_w_border <- as.factor(nn.model.data$interface_w_border)
nn.model.data.h2o <- h2o::as.h2o(nn.model.data[,-c(1)])

#split into train and validate 
data.splits <- h2o::h2o.splitFrame(data =  nn.model.data.h2o, ratios = 0.8, seed = 1234)
train <- data.splits[[1]]
valid <- data.splits[[2]]

#train model 
fit.dl = h2o::h2o.deeplearning( x = names(train)[-which(names(train) == 'interface_w_border')], y =  names(train)[which(names(train) == 'interface_w_border')], 
                                training_frame = train, validation_frame = valid, 
                                activation = "Tanh", autoencoder = F,
                                hidden = c(400,200,100,10), epochs = 1000,
                                export_weights_and_biases = TRUE )

#Get performance on training data 
train.performance <- h2o::h2o.performance(fit.dl, train = T)
sink(paste0(out.dir, "/training_results.txt"))
print(train.performance)
sink()

#Extract confusion matrix
cm_raw <- as.data.frame(h2o::h2o.confusionMatrix(train.performance))
cm_matrix <- cm_raw[1:3, 1:3]
rownames(cm_matrix) <- cm_raw[1:3, "Row.labels"]
cm_melted <- melt(as.matrix(cm_matrix))
colnames(cm_melted) <- c("Actual", "Predicted", "Count")

#Generate heatmap
ggplot(cm_melted, aes(x = Predicted, y = Actual, fill = Count)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Count), color = "black", size = 5, fontface = "bold") +
  scale_fill_gradient(low = "#e6f2ff", high = "#0066cc") +
  scale_y_discrete(limits = c("Active", "Intermediate", "Passive"), 
                   labels = c("Active", "Intermediate", "Passive")) +
  labs(title = "H2O Deep Learning Confusion Matrix",
       x = "Predicted Class",
       y = "Actual Class") +
  theme_minimal(base_size = 14) +
  theme(panel.grid = element_blank())
ggsave(paste0(out.dir, '/Confusion_matrix.png'))

#Plot feature importances
png(paste0(out.dir, "variable_importance.png"), width = 1000, height = 600, res = 120)
h2o::h2o.varimp_plot(fit.dl, num_of_features = 10)
dev.off()

#create prediction object 
prediction.data <- input.data |> dplyr::filter(interface_w_border == 'Interface')
colnames(prediction.data)[1] <- 'Barcodes'
prediction.h2o <- h2o::as.h2o(prediction.data[,-which(names(prediction.data) %in% c("Barcodes", "interface_w_border"))])

#Get prediction 
prediction.res <- h2o::h2o.predict(fit.dl, prediction.h2o)

#Add new annotations to original data 
input.data <- as.data.frame(input.data, row.names = input.data$Barcodes)
colnames(input.data)[1] <- 'Barcodes'
input.data$model.annotations <- input.data$interface_w_border
input.data[input.data$Barcodes %in% prediction.data$Barcodes,]$model.annotations <- as.vector(prediction.res[,1])

#shutdown conection
h2o::h2o.shutdown(prompt = F) 
if("package:h2o" %in% search()) {
  detach("package:h2o", unload = TRUE)
}


#save annotations 
object.subset$model_predictions_interface <- input.data$model.annotations
interface.anno <- as.data.frame(object.subset$model_predictions_interface)
colnames(interface.anno) <- "global_tri_state"
write.csv(interface.anno, paste0(out.dir, '/interface_predictions_annotations.csv'))

#load data
interface.anno <- read.csv(paste0(out.dir, '/interface_predictions_annotations.csv'))
object.interface$model_predictions_interface <- interface.anno$global_tri_state


#save plots and annotation
#SpatialPlot(object.subset, 'model_predictions_interface',pt.size.factor = size[sample.num],
#                   image.alpha = 0.3, cols = cols)
#ggsave(paste0(out.dir, '/Interface_annotations.png'))


### Make plots for interface
for(sample.num in 1:9){
  #select sample
  sample <- letters[sample.num]
  
  #make output dir
  output.dir <- paste0(out.dir, sample)
  dir.create(output.dir)
  
  #subset samples
  object.sample <- subset(object.subset , subset = orig.ident == paste0('sample_', sample))
  
  #Plots predicted interface
  SpatialPlot(object.sample , 'model_predictions_interface',pt.size.factor = size[sample.num],
              image.alpha = 0.3, cols = cols)
  ggsave(paste0(output.dir, '/cluster_activity_plot.png'))
  
  VlnPlot(object.sample, 'Active_Score1', group.by = 'model_predictions_interface', cols = cols)
  ggsave(paste0(output.dir, '/activity_vln_plot.png'))
  
  #save annotation
  write.csv(as.data.frame(object.sample$model_predictions_interface), paste0(output.dir, '/predicted_annotations.csv'),
            quote = F)
  
}

# Read back annotation on sample 
object.interface.sub <- subset(object.interface, subset = orig.ident %in% paste0('sample_', c('A', 'E', 'G', 'I', 'K', 'M')))
anno <- c()
for(i in  c('A', 'E', 'G', 'I', 'K', 'M')){
  sample.anno <- read.csv(paste0(base.dir, sample, '/interface_res/predicted_annotations.csv'))
  anno <- c(anno, sample.anno[,2])
}
object.interface.sub$interface.annotation <- anno


##### 2. Plot Variogram of Interface By Sample ########################################

#subset by sample 
coords.list <- list()
scores.list <- list()
for(i in letters){
  object.sample <- subset(object.interface, subset = orig.ident == paste0('sample_', i))
  coords.list[[i]] <- GetTissueCoordinates(object.sample)
  scores.list[[i]] <-  object.sample$Active_Score1
}

#Make variogram
par(mar = c(5, 5, 4, 2) + 0.1, cex.main = 1, cex.lab = 0.9)
for(i in letters){
  #output.dir <- paste0(out.dir, i)
  sp.tme=coords.list[[i]]
  #cell.tme=as.factor(scores.list[[i]])
  cell.tme=scores.list[[i]]

  #combine df
  comb <- data.frame(
    imagerow = sp.tme$imagerow,
    imagecol = sp.tme$imagecol,
    score = cell.tme
  )
  comb <- na.omit(comb) 
  
  #Get scale for micrometers
  pixel_dist_matrix <- as.matrix(dist(comb[, c("imagerow", "imagecol")]))
  diag(pixel_dist_matrix) <- Inf
  microns_per_pixel <- 100 / min(pixel_dist_matrix)
  
  #Convert coordinates to micrometers
  comb$x_um <- comb$imagecol * microns_per_pixel
  comb$y_um <- comb$imagerow * microns_per_pixel
  
  #Compute variogram
  coordinates(comb) <- ~ x_um + y_um
  v_obj <- variogram(score ~ 1, comb)
  
  #Plot
  #png(paste0(output.dir, "/activity_variogram.png"), width = 600, height = 600, res = 120)
  plot(v_obj$dist, v_obj$gamma, 
       main = paste0("Activity Score: Sample", i), 
       xlab = expression("distance (" * mu * "m)"), 
       ylab = "semivariance", 
       pch=19, cex.lab=2, cex.axis=2, cex=2, cex.main=2)
  #dev.off()
}


##### 3. Calculate autocorrelation of spots in Interface Region ####################
scores.list <- list()
for(i in letters){
  object.sample <- subset(object.subset, subset = orig.ident == paste0('sample_', i))
  factor.mat <- data.frame(category = as.factor(object.sample$model_predictions_interface))
  one.hot.mat <- model.matrix(~ category - 1, data = factor.mat)
  colnames(one.hot.mat) <- gsub('category', '', colnames(one.hot.mat))
  scores.list[[i]] <-  data.frame(Active_score = object.sample$Active_Score1, one.hot.mat)
}
autocor.res <- scores.autocor(coords.list, scores.list)

autocor.res.p.values <- scores.autocor(coords.list, scores.list, return.pvalue = T)

#Plot Heatmap
Heatmap(autocor.res, name = "Moran's I", column_title = 'Autocorrelation at Extended Interface')

#save results 
write.csv(autocor.res, paste0(base.dir, 'interface_res/autcor_results.csv'))
write.csv(autocor.res.p.values, paste0(base.dir, 'interface_res/autcor_results_pvalues.csv'))


##### 4. Get the Area of activity niches  #######################################

#etract active niche object
active.niche.object <- subset(object.interface, subset = model_predictions_interface == 'Active')

#run through every sample and calculate size
area.list <- list()
for(sample in letters){
  active.niche.object.subset <- subset(active.niche.object, subset = orig.ident == paste0('sample_', sample))
  
  #extract coordinates and annotation
  coords <- GetTissueCoordinates(active.niche.object.subset)
  
  #Create graph
  sp.tme <- coords 
  colnames(sp.tme) <- c('row', 'col')
  pt.dist <- min(dist(sp.tme))
  sp::coordinates(sp.tme) <- ~ row  + col
  neib <- spdep::dnearneigh(as.matrix(sp::coordinates(sp.tme)), d1 = 1, d2 = pt.dist + pt.dist/2)
  
  #Remove niehgborhoods with less than 3 spots 
  valid_pts <- card(neib) >=1
  sp.tme_cleaned <- sp.tme[valid_pts, ]
  neib <- dnearneigh(
    as.matrix(sp::coordinates(sp.tme_cleaned)), 
    d1 = 1, 
    d2 = pt.dist + pt.dist/2
  )
  
  #Clean coordinates
  coords_clean <- as.matrix(sp::coordinates(sp.tme_cleaned))
  
  #Get areas
  areas <- sapply(seq_along(neib), function(i) {
    neighbor_indices <- c(i, neib[[i]])
    
    #Extract points and compute convex hull
    pts <- coords_clean[neighbor_indices, , drop = FALSE]
    hull_idx <- chull(pts)
    hull_pts <- pts[c(hull_idx, hull_idx[1]), ]
    poly <- st_polygon(list(as.matrix(hull_pts)))
    
    return(st_area(poly))
  })
  
  #add to list
  area.list[[sample]] <- areas
}


#Subset for sample G
sample <- 'G'
area.vec <- unlist(area.list$A) 
active.niche.object.subset <- subset(active.niche.object, subset = orig.ident == paste0('sample_', sample))
coords <- Seurat::GetTissueCoordinates(active.niche.object.subset, scale.factor = "lowres")

#Visulize area 
SpatialPlot(active.niche.object.subset, group.by = 'model_predictions_interface', pt.size.factor = size[sample] + 100)

#Plot density of areas
df_areas <- bind_rows(
  lapply(names(areas.list), function(sample_id) {
    data.frame(
      sample = sample_id,
      area = unname(areas.list[[sample_id]])
    )
  })
)


ggplot(df_areas, aes(x = area, fill = sample, color = sample)) +
  geom_density(alpha = 0.2, linewidth = 0.8) +
  scale_x_log10(labels = scales::comma) +
  theme_classic() +
  labs(
    title = "Active Niche Area Density Distribution Across Samples",
    x = expression("Area (" * mu * "m"^2 * ") - Log10 Scale"),
    y = "Density",
    fill = "Sample",
    color = "Sample"
  )

#Get mean of areas
mean.list <- c()
for(sample in letters){
  mean.list <- c(mean.list, mean(areas.list[[sample]]))
}


