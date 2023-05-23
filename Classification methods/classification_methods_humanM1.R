<<<<<<< HEAD
################################
##### Seurat sample usage ######
################################

library(Seurat)
library(Matrix)
library(ggplot2)
library(knitr)
library(dplyr)
set.seed(3456)
file_path_root = "/root/datasets/humanM1/"
train_data = read.csv(paste0(file_path_root, "data1/data.csv"),row.names=1)
train_meta = read.csv(paste0(file_path_root, "data1/meta.csv"),row.names=1)
colnames(train_data) = rownames(train_meta)
test_data = read.csv(paste0(file_path_root, "data2/data.csv"),row.names=1)
test_meta = read.csv(paste0(file_path_root, "data2/meta.csv"),row.names=1)
colnames(test_data) = rownames(test_meta)
#load reference
reference_count <- train_data
reference_meta <- train_meta
reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "example")
reference_seurat <- AddMetaData(reference_seurat, reference_meta)

#load query
query_count <- test_data
query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "example")

#standard pipeline
reference_seurat <- NormalizeData(object = reference_seurat)
reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = 'vst', nfeatures = 2000)
reference_seurat <- ScaleData(reference_seurat)

query_seurat <- NormalizeData(object = query_seurat)
query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = 'vst', nfeatures = 2000)
query_seurat <- ScaleData(query_seurat)

##prediction###
sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = query_seurat,
                                   dims = 1:30)
##replace Group with the actual column name from meta
predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$subclass_label,
                            dims = 1:30)
query_seurat <- AddMetaData(object = query_seurat, metadata = predictions)

predict_meta = data.frame(rownames(query_seurat@meta.data), query_seurat@meta.data$predicted.id)

rownames(predict_meta) = predict_meta$rownames.query_seurat.meta.data.
predict_meta = predict_meta[rownames(test_meta),]
true_label = test_meta$subclass_label
sum(predict_meta$query_seurat.meta.data.predicted.id==true_label) / length(true_label)


saveRDS(query_seurat, paste0(file_path_root,"seurat_result.rds"))


################################
##### SingleR sample usage ######
################################

library(Matrix)
library(SingleR)
library(ggplot2)
library(knitr)
library(genefilter)
set.seed(3456)




predictions <- SingleR(test=test_data, assay.type.test=1, 
                       ref=train_data, labels=train_meta$subclass_label)


table(predictions$labels)
#sum(test_meta$subclass_label==predictions$labels)/length(predictions$labels)

saveRDS(predictions, paste0(file_path_root,"singleR_result.rds"))
#######################################
##### scmap sample usage ######
#######################################
library(SingleCellExperiment)
library(scmap)
library(Matrix)
set.seed(3456)


# Load reference meta
reference_meta <- train_meta
# load reference data
reference_norm <- train_data
# replace Group with the right column name in meta data
ref_ann <-  as.data.frame(reference_meta$subclass_label)
colnames(ref_ann) <- "celltype"
# Create sce object
reference <- SingleCellExperiment(assays = list(normcounts = as.matrix(reference_norm)), colData = ref_ann)
logcounts(reference) <- normcounts(reference)
rowData(reference)$feature_symbol <- rownames(reference)
reference <- reference[!duplicated(rownames(reference)), ]
reference <- selectFeatures(reference, suppress_plot = FALSE)

#scmap-cell
reference <- indexCell(reference)

#query data
query_norm <- test_data
query <- SingleCellExperiment(assays = list(normcounts = as.matrix(query_norm)))
logcounts(query) <- normcounts(query)
rowData(query)$feature_symbol <- rownames(query)

# Prediction
scmapCell_results <- scmapCell(
  projection = query,
  list(
    ref = metadata(reference)$scmap_cell_index
  )
)

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(reference)$celltype)
  ), threshold = -Inf
)

result <- scmapCell_clusters$scmap_cluster_labs
row.names(result) <-  colnames(scmapCell_results$ref$cells)

#sum(test_meta$subclass_label==result[,1])/length(result[,1])


saveRDS(result, paste0(file_path_root,"scmap_result.rds"))








#######################################
##### singleCellNet sample usage ######
#######################################
library(singleCellNet)
library(dplyr)
set.seed(3456)




###For smaller dataset: cell type smaller than 4 cannot be split in singleCellNet and need to be exclude out#####
splitCommon<-function(sampTab, ncells, dLevel="cell_ontology_class"){
  cts<-unique(as.vector(sampTab[,dLevel]))
  trainingids<-vector()
  for(ct in cts){
    stX<-sampTab[sampTab[,dLevel]==ct, ,drop = F]
    if(nrow(stX) > 4){
      cat(ct, ": ")
      ccount<-nrow(stX)-3
      ccount<-min(ccount, ncells)
      cat(nrow(stX),"\n")
      trainingids<-append(trainingids, sample(rownames(stX), ccount))
    }
  }
  val_ids<-setdiff(rownames(sampTab), trainingids)
  list(train=sampTab[trainingids, , drop = F], val=sampTab[val_ids, , drop = F])
}
####################################################################################################################

#load query
query_count <- test_data
query_meta <- test_meta
gene_qeury <- rownames(query_count)

#load referecne/train
reference_count <- train_data
reference_meta <- train_meta
reference_meta <- droplevels(reference_meta)

#process ref
commonGenes <- intersect(rownames(reference_count), gene_qeury)
reference_count <- reference_count[commonGenes, ]

#Split for training and assessment, and transform training data
#dLevel need to set to the correct column name for meta data
reference_list <- splitCommon(reference_meta, ncells = 500, dLevel = "subclass_label")
reference_train <- reference_list[[1]]
expTrain <- reference_count[, rownames(reference_train)]
expTrain = as.matrix(expTrain)





######################
class_info<-scn_train(stTrain = reference_train, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "subclass_label")

classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=query_count, nrand = 50)
result <- c()
for(i in 1:ncol(classRes_val_all)){
  result <- c(result,rownames(classRes_val_all)[(classRes_val_all[,i]==max(classRes_val_all[,i]))][1])
}
names(result) = colnames(classRes_val_all)
result = result[rownames(test_meta)]

#sum(test_meta$subclass_label==result)/length(result)
saveRDS(result, paste0(file_path_root,"singleCellNet_result.rds"))
#saveRDS(result, "/picb/neurosys/chiyuhao/PFC/final_human_PFC/singleCellNet_result.rds")




################################
##### CHETAH sample usage ######
################################
library(CHETAH)
set.seed(3456)


# Load meta
ref_meta <- train_meta

## Make SingleCellExperiments
### reference need to be normalized
reference_ct <- train_data
### replace Group with correct column name in meta
ref_ct <- as.character(ref_meta$subclass_label)
reference <- SingleCellExperiment(assays = list(counts = reference_ct),
                                  colData = DataFrame(celltypes = ref_ct))

###optional normalized for query
query_ct <- test_data
input <- SingleCellExperiment(assays = list(counts = query_ct))

## Run CHETAH, minimize unassign by set threshold as -Inf by thresh = -Inf, otherwise use default parameter
input <- CHETAHclassifier(input = input, ref_cells = reference)
## Extract celltypes
CHETAH_pred <- as.data.frame(input$celltype_CHETAH)
colnames(CHETAH_pred) <- 'CHETAH_pred'
CHETAH_pred = CHETAH_pred[rownames(test_meta),]
#sum(test_meta$subclass_label==CHETAH_pred)/length(CHETAH_pred)
#saveRDS(CHETAH_pred, "/picb/neurosys/chiyuhao/PFC/final_human_PFC/CHETAH_result.rds")
=======
################################
##### Seurat sample usage ######
################################

library(Seurat)
library(Matrix)
library(ggplot2)
library(knitr)
library(dplyr)
set.seed(3456)
file_path_root = "/root/datasets/humanM1/"
train_data = read.csv(paste0(file_path_root, "data1/data.csv"),row.names=1)
train_meta = read.csv(paste0(file_path_root, "data1/meta.csv"),row.names=1)
colnames(train_data) = rownames(train_meta)
test_data = read.csv(paste0(file_path_root, "data2/data.csv"),row.names=1)
test_meta = read.csv(paste0(file_path_root, "data2/meta.csv"),row.names=1)
colnames(test_data) = rownames(test_meta)
#load reference
reference_count <- train_data
reference_meta <- train_meta
reference_seurat <- CreateSeuratObject(counts = reference_count, min.cells = 0, min.features = 0, project = "example")
reference_seurat <- AddMetaData(reference_seurat, reference_meta)

#load query
query_count <- test_data
query_seurat <- CreateSeuratObject(counts = query_count, min.cells = 0, min.features = 0, project = "example")

#standard pipeline
reference_seurat <- NormalizeData(object = reference_seurat)
reference_seurat <- FindVariableFeatures(object = reference_seurat, selection.method = 'vst', nfeatures = 2000)
reference_seurat <- ScaleData(reference_seurat)

query_seurat <- NormalizeData(object = query_seurat)
query_seurat <- FindVariableFeatures(object = query_seurat, selection.method = 'vst', nfeatures = 2000)
query_seurat <- ScaleData(query_seurat)

##prediction###
sim.anchors <- FindTransferAnchors(reference = reference_seurat, query = query_seurat,
                                   dims = 1:30)
##replace Group with the actual column name from meta
predictions <- TransferData(anchorset = sim.anchors, refdata = reference_seurat$subclass_label,
                            dims = 1:30)
query_seurat <- AddMetaData(object = query_seurat, metadata = predictions)

predict_meta = data.frame(rownames(query_seurat@meta.data), query_seurat@meta.data$predicted.id)

rownames(predict_meta) = predict_meta$rownames.query_seurat.meta.data.
predict_meta = predict_meta[rownames(test_meta),]
true_label = test_meta$subclass_label
sum(predict_meta$query_seurat.meta.data.predicted.id==true_label) / length(true_label)


saveRDS(query_seurat, paste0(file_path_root,"seurat_result.rds"))


################################
##### SingleR sample usage ######
################################

library(Matrix)
library(SingleR)
library(ggplot2)
library(knitr)
library(genefilter)
set.seed(3456)




predictions <- SingleR(test=test_data, assay.type.test=1, 
                       ref=train_data, labels=train_meta$subclass_label)


table(predictions$labels)
#sum(test_meta$subclass_label==predictions$labels)/length(predictions$labels)

saveRDS(predictions, paste0(file_path_root,"singleR_result.rds"))
#######################################
##### scmap sample usage ######
#######################################
library(SingleCellExperiment)
library(scmap)
library(Matrix)
set.seed(3456)


# Load reference meta
reference_meta <- train_meta
# load reference data
reference_norm <- train_data
# replace Group with the right column name in meta data
ref_ann <-  as.data.frame(reference_meta$subclass_label)
colnames(ref_ann) <- "celltype"
# Create sce object
reference <- SingleCellExperiment(assays = list(normcounts = as.matrix(reference_norm)), colData = ref_ann)
logcounts(reference) <- normcounts(reference)
rowData(reference)$feature_symbol <- rownames(reference)
reference <- reference[!duplicated(rownames(reference)), ]
reference <- selectFeatures(reference, suppress_plot = FALSE)

#scmap-cell
reference <- indexCell(reference)

#query data
query_norm <- test_data
query <- SingleCellExperiment(assays = list(normcounts = as.matrix(query_norm)))
logcounts(query) <- normcounts(query)
rowData(query)$feature_symbol <- rownames(query)

# Prediction
scmapCell_results <- scmapCell(
  projection = query,
  list(
    ref = metadata(reference)$scmap_cell_index
  )
)

scmapCell_clusters <- scmapCell2Cluster(
  scmapCell_results, 
  list(
    as.character(colData(reference)$celltype)
  ), threshold = -Inf
)

result <- scmapCell_clusters$scmap_cluster_labs
row.names(result) <-  colnames(scmapCell_results$ref$cells)

#sum(test_meta$subclass_label==result[,1])/length(result[,1])


saveRDS(result, paste0(file_path_root,"scmap_result.rds"))








#######################################
##### singleCellNet sample usage ######
#######################################
library(singleCellNet)
library(dplyr)
set.seed(3456)




###For smaller dataset: cell type smaller than 4 cannot be split in singleCellNet and need to be exclude out#####
splitCommon<-function(sampTab, ncells, dLevel="cell_ontology_class"){
  cts<-unique(as.vector(sampTab[,dLevel]))
  trainingids<-vector()
  for(ct in cts){
    stX<-sampTab[sampTab[,dLevel]==ct, ,drop = F]
    if(nrow(stX) > 4){
      cat(ct, ": ")
      ccount<-nrow(stX)-3
      ccount<-min(ccount, ncells)
      cat(nrow(stX),"\n")
      trainingids<-append(trainingids, sample(rownames(stX), ccount))
    }
  }
  val_ids<-setdiff(rownames(sampTab), trainingids)
  list(train=sampTab[trainingids, , drop = F], val=sampTab[val_ids, , drop = F])
}
####################################################################################################################

#load query
query_count <- test_data
query_meta <- test_meta
gene_qeury <- rownames(query_count)

#load referecne/train
reference_count <- train_data
reference_meta <- train_meta
reference_meta <- droplevels(reference_meta)

#process ref
commonGenes <- intersect(rownames(reference_count), gene_qeury)
reference_count <- reference_count[commonGenes, ]

#Split for training and assessment, and transform training data
#dLevel need to set to the correct column name for meta data
reference_list <- splitCommon(reference_meta, ncells = 500, dLevel = "subclass_label")
reference_train <- reference_list[[1]]
expTrain <- reference_count[, rownames(reference_train)]
expTrain = as.matrix(expTrain)





######################
class_info<-scn_train(stTrain = reference_train, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "subclass_label")

classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=query_count, nrand = 50)
result <- c()
for(i in 1:ncol(classRes_val_all)){
  result <- c(result,rownames(classRes_val_all)[(classRes_val_all[,i]==max(classRes_val_all[,i]))][1])
}
names(result) = colnames(classRes_val_all)
result = result[rownames(test_meta)]

#sum(test_meta$subclass_label==result)/length(result)
saveRDS(result, paste0(file_path_root,"singleCellNet_result.rds"))
#saveRDS(result, "/picb/neurosys/chiyuhao/PFC/final_human_PFC/singleCellNet_result.rds")




################################
##### CHETAH sample usage ######
################################
library(CHETAH)
set.seed(3456)


# Load meta
ref_meta <- train_meta

## Make SingleCellExperiments
### reference need to be normalized
reference_ct <- train_data
### replace Group with correct column name in meta
ref_ct <- as.character(ref_meta$subclass_label)
reference <- SingleCellExperiment(assays = list(counts = reference_ct),
                                  colData = DataFrame(celltypes = ref_ct))

###optional normalized for query
query_ct <- test_data
input <- SingleCellExperiment(assays = list(counts = query_ct))

## Run CHETAH, minimize unassign by set threshold as -Inf by thresh = -Inf, otherwise use default parameter
input <- CHETAHclassifier(input = input, ref_cells = reference)
## Extract celltypes
CHETAH_pred <- as.data.frame(input$celltype_CHETAH)
colnames(CHETAH_pred) <- 'CHETAH_pred'
CHETAH_pred = CHETAH_pred[rownames(test_meta),]
#sum(test_meta$subclass_label==CHETAH_pred)/length(CHETAH_pred)
#saveRDS(CHETAH_pred, "/picb/neurosys/chiyuhao/PFC/final_human_PFC/CHETAH_result.rds")
>>>>>>> 3e2f6a6 (add large file)
saveRDS(CHETAH_pred, paste0(file_path_root,"CHETAH_result.rds"))