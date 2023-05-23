<<<<<<< HEAD
#consensus 1
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
library(CYANO)
file_path_root = "/picb/neurosys/chiyuhao/PFC/mouse_M1/34616066/10X_cells_v2_AIBS/data2/"
de.param <- de_param(padj.th=0.01,
                     q1.th       = 0.4, 
                     q.diff.th   = 0.7, 
                     de.score.th = 150,
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
cell_number=10
new_data  = data
result_1 = readRDS(paste0(file_path_root, "consesus_1_result_1.rds"))
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
ref_meta = read.csv(paste0(file_path_root, "meta.csv"),row.names=1)
rownames(ref_meta) <- gsub("-","\\.",rownames(ref_meta))
anno_meta = RunSubclassClassify(new_data, group_meta, new_data, ref_meta)
saveRDS(anno_meta, paste0(file_path_root, "classificiation_resault.rds"))
new_data = new_data[,rownames(anno_meta)]
result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",TRUE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "pse_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "pse_roc_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "pse_roc_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_roc.csv"),
=======
#consensus 1
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
library(CYANO)
file_path_root = "/picb/neurosys/chiyuhao/PFC/mouse_M1/34616066/10X_cells_v2_AIBS/data2/"
de.param <- de_param(padj.th=0.01,
                     q1.th       = 0.4, 
                     q.diff.th   = 0.7, 
                     de.score.th = 150,
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
cell_number=10
new_data  = data
result_1 = readRDS(paste0(file_path_root, "consesus_1_result_1.rds"))
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
ref_meta = read.csv(paste0(file_path_root, "meta.csv"),row.names=1)
rownames(ref_meta) <- gsub("-","\\.",rownames(ref_meta))
anno_meta = RunSubclassClassify(new_data, group_meta, new_data, ref_meta)
saveRDS(anno_meta, paste0(file_path_root, "classificiation_resault.rds"))
new_data = new_data[,rownames(anno_meta)]
result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",TRUE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "pse_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "pse_roc_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "pse_roc_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_roc.csv"),
>>>>>>> 3e2f6a6 (add large file)
                   paste0(file_path_root, "pse_roc_old.csv"),paste0(file_path_root, "pse_roc_raw.csv"))