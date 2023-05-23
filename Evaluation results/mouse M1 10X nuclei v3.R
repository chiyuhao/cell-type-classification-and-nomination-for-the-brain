<<<<<<< HEAD
library(scrattch.hicat)
library(Matrix)
library(matrixStats)
library(CYANO)
library(ROGUE)
library(tibble)
library(ktplots)
library(beeswarm)
library(WGCNA)
library(edgeR)   
library(feather)
library(dendextend)
library(ggplot2)
library(dplyr)
library(Seurat)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
library(knitr)
library(parallel)
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)


source("/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/get_file_name_function_old.R")
file_path_root = "/picb/neurosys/chiyuhao/PFC/mouse_M1/34616066/10X_nuclei_v2_AIBS/"
de.param <- de_param(padj.th=0.01,
                     q1.th       = 0.4, 
                     q.diff.th   = 0.7, 
                     de.score.th = 100,
                     min.cells=30)


data = read.csv(paste0(file_path_root,"data.csv"),row.names=1)
start <- Sys.time()
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"run_cluster"))
result_1 <- RunConsensus1(norm.dat,
                          de.param = de.param,
                          override = TRUE,
                          output_dir = paste0(file_path_root,"run_cluster"),
                          mc.cores = 1)
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
saveRDS(result_1, paste0(file_path_root,"consesus_1_result_1.rds")) 
saveRDS(group_meta, paste0(file_path_root,"group_meta_1_new.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root,"running_time_consensus1.rds"))



new_sample_data=readRDS(paste0(file_path_root,"pse_raw_data_1.rds"))
rownames(new_sample_data) <- 1:nrow(new_sample_data)
new_sample_data = t(new_sample_data)
new_meta = read.csv(paste0(file_path_root, "new_meta_pse_1.csv"),row.names=1)
mouse_data <- CreateSeuratObject(counts = new_sample_data, min.cells = 0, min.features = 0, project = "example")
mouse_data <- AddMetaData(mouse_data, new_meta)
mouse_data <- NormalizeData(mouse_data, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(mouse_data) <- mouse_data$subclass_label
all_marker_list = list()
num = 1
for(i in unique(mouse_data$subclass_label)){
  temp_mouse_data = subset(mouse_data, idents = i)
  Idents(temp_mouse_data) <- temp_mouse_data$cluster_label
  #plan(workers = 6)
  mouse_cells_markers <- FindAllMarkers(temp_mouse_data, test.use = "roc",densify=T)
  mouse_cells_markers = mouse_cells_markers[mouse_cells_markers$avg_log2FC>0,]
  all_marker_list[[num]] = mouse_cells_markers
  num = num + 1
}
names(all_marker_list) = unique(mouse_data$subclass_label)
pse_meta = new_meta
pse_meta$class_label = pse_meta$subclass_label
pse_data = new_sample_data
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_roc.csv"),
                   paste0(file_path_root, "pse_roc_old.csv"),paste0(file_path_root, "pse_roc_raw.csv"))
=======
library(scrattch.hicat)
library(Matrix)
library(matrixStats)
library(CYANO)
library(ROGUE)
library(tibble)
library(ktplots)
library(beeswarm)
library(WGCNA)
library(edgeR)   
library(feather)
library(dendextend)
library(ggplot2)
library(dplyr)
library(Seurat)
library(matrixStats)
library(Matrix)
library(scrattch.hicat)
library(knitr)
library(parallel)
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)


source("/picb/neurosys/chiyuhao/PFC/final_datasets/datasets/get_file_name_function_old.R")
file_path_root = "/picb/neurosys/chiyuhao/PFC/mouse_M1/34616066/10X_nuclei_v2_AIBS/"
de.param <- de_param(padj.th=0.01,
                     q1.th       = 0.4, 
                     q.diff.th   = 0.7, 
                     de.score.th = 100,
                     min.cells=30)


data = read.csv(paste0(file_path_root,"data.csv"),row.names=1)
start <- Sys.time()
data = as.matrix(data)
norm.dat <- Matrix(cpm(data), sparse = TRUE)
norm.dat <- log2(cpm(norm.dat)+1)
set.seed(3456)


dir.create(paste0(file_path_root,"run_cluster"))
result_1 <- RunConsensus1(norm.dat,
                          de.param = de.param,
                          override = TRUE,
                          output_dir = paste0(file_path_root,"run_cluster"),
                          mc.cores = 1)
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
saveRDS(result_1, paste0(file_path_root,"consesus_1_result_1.rds")) 
saveRDS(group_meta, paste0(file_path_root,"group_meta_1_new.rds")) 
end <- Sys.time()
runningtime <- end-start
saveRDS(runningtime, paste0(file_path_root,"running_time_consensus1.rds"))



new_sample_data=readRDS(paste0(file_path_root,"pse_raw_data_1.rds"))
rownames(new_sample_data) <- 1:nrow(new_sample_data)
new_sample_data = t(new_sample_data)
new_meta = read.csv(paste0(file_path_root, "new_meta_pse_1.csv"),row.names=1)
mouse_data <- CreateSeuratObject(counts = new_sample_data, min.cells = 0, min.features = 0, project = "example")
mouse_data <- AddMetaData(mouse_data, new_meta)
mouse_data <- NormalizeData(mouse_data, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(mouse_data) <- mouse_data$subclass_label
all_marker_list = list()
num = 1
for(i in unique(mouse_data$subclass_label)){
  temp_mouse_data = subset(mouse_data, idents = i)
  Idents(temp_mouse_data) <- temp_mouse_data$cluster_label
  #plan(workers = 6)
  mouse_cells_markers <- FindAllMarkers(temp_mouse_data, test.use = "roc",densify=T)
  mouse_cells_markers = mouse_cells_markers[mouse_cells_markers$avg_log2FC>0,]
  all_marker_list[[num]] = mouse_cells_markers
  num = num + 1
}
names(all_marker_list) = unique(mouse_data$subclass_label)
pse_meta = new_meta
pse_meta$class_label = pse_meta$subclass_label
pse_data = new_sample_data
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_roc.csv"),
                   paste0(file_path_root, "pse_roc_old.csv"),paste0(file_path_root, "pse_roc_raw.csv"))
>>>>>>> 3e2f6a6 (add large file)
