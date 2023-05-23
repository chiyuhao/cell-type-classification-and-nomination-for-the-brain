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


source("/root/datasets/get_file_name_function.R")

cell_number=2
file_path_root = "/root/datasets/humanPFC/data1/"
new_data  = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
result_1 = readRDS(paste0(file_path_root, "consesus_1_result_1.rds"))
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
ref_meta = read.csv(paste0(file_path_root, "meta.csv"),row.names=1)
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
                   paste0(file_path_root, "pse_roc_old.csv"),paste0(file_path_root, "pse_roc_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",TRUE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "pse_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "pse_wilcox_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "pse_wilcox_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_wilcox.csv"),
                   paste0(file_path_root, "pse_wilcox_old.csv"),paste0(file_path_root, "pse_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_roc_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_roc.csv"),
                   paste0(file_path_root, "all_roc_old.csv"),paste0(file_path_root, "all_roc_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_wilcox_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_wilcox.csv"),
                   paste0(file_path_root, "all_wilcox_old.csv"),paste0(file_path_root, "all_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"t",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_t.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_t_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_t.csv"),
                   paste0(file_path_root, "all_t_old.csv"),paste0(file_path_root, "all_t_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"LR",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_lr.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_lr_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_lr.csv"),
                   paste0(file_path_root, "all_lr_old.csv"),paste0(file_path_root, "all_lr_raw.csv"))





file_path_root = "/root/datasets/humanPFC/data2/"
new_data  = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
result_1 = readRDS(paste0(file_path_root, "consesus_1_result_2.rds"))
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
ref_meta = read.csv(paste0(file_path_root, "meta.csv"),row.names=1)
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
                   paste0(file_path_root, "pse_roc_old.csv"),paste0(file_path_root, "pse_roc_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",TRUE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "pse_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "pse_wilcox_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "pse_wilcox_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_wilcox.csv"),
                   paste0(file_path_root, "pse_wilcox_old.csv"),paste0(file_path_root, "pse_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_roc_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_roc.csv"),
                   paste0(file_path_root, "all_roc_old.csv"),paste0(file_path_root, "all_roc_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_wilcox_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_wilcox.csv"),
                   paste0(file_path_root, "all_wilcox_old.csv"),paste0(file_path_root, "all_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"t",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_t.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_t_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_t.csv"),
                   paste0(file_path_root, "all_t_old.csv"),paste0(file_path_root, "all_t_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"LR",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_lr.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_lr_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_lr.csv"),
                   paste0(file_path_root, "all_lr_old.csv"),paste0(file_path_root, "all_lr_raw.csv"))
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


source("/root/datasets/get_file_name_function.R")

cell_number=2
file_path_root = "/root/datasets/humanPFC/data1/"
new_data  = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
result_1 = readRDS(paste0(file_path_root, "consesus_1_result_1.rds"))
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
ref_meta = read.csv(paste0(file_path_root, "meta.csv"),row.names=1)
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
                   paste0(file_path_root, "pse_roc_old.csv"),paste0(file_path_root, "pse_roc_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",TRUE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "pse_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "pse_wilcox_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "pse_wilcox_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_wilcox.csv"),
                   paste0(file_path_root, "pse_wilcox_old.csv"),paste0(file_path_root, "pse_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_roc_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_roc.csv"),
                   paste0(file_path_root, "all_roc_old.csv"),paste0(file_path_root, "all_roc_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_wilcox_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_wilcox.csv"),
                   paste0(file_path_root, "all_wilcox_old.csv"),paste0(file_path_root, "all_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"t",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_t.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_t_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_t.csv"),
                   paste0(file_path_root, "all_t_old.csv"),paste0(file_path_root, "all_t_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"LR",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_lr.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_lr_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_lr.csv"),
                   paste0(file_path_root, "all_lr_old.csv"),paste0(file_path_root, "all_lr_raw.csv"))





file_path_root = "/root/datasets/humanPFC/data2/"
new_data  = read.csv(paste0(file_path_root, "data.csv"),row.names=1)
result_1 = readRDS(paste0(file_path_root, "consesus_1_result_2.rds"))
group_meta = data.frame(names(result_1$cl),result_1$cl)
colnames(group_meta) = c("cell","group")
ref_meta = read.csv(paste0(file_path_root, "meta.csv"),row.names=1)
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
                   paste0(file_path_root, "pse_roc_old.csv"),paste0(file_path_root, "pse_roc_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",TRUE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "pse_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "pse_wilcox_meta.rds"))
saveRDS(pse_data, paste0(file_path_root, "pse_wilcox_data.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "pse_wilcox.csv"),
                   paste0(file_path_root, "pse_wilcox_old.csv"),paste0(file_path_root, "pse_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"roc",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_roc.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_roc_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_roc.csv"),
                   paste0(file_path_root, "all_roc_old.csv"),paste0(file_path_root, "all_roc_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"wilcox",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_wilcox.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_wilcox_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_wilcox.csv"),
                   paste0(file_path_root, "all_wilcox_old.csv"),paste0(file_path_root, "all_wilcox_raw.csv"))


result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"t",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_t.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_t_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_t.csv"),
                   paste0(file_path_root, "all_t_old.csv"),paste0(file_path_root, "all_t_raw.csv"))

result = RunFindDEGeneNew(new_data, anno_meta,cell_number,"LR",FALSE)
all_marker_list = result[[1]]
saveRDS(all_marker_list, paste0(file_path_root, "all_lr.rds"))
pse_meta = result[[2]]
pse_data = result[[3]]
saveRDS(pse_meta, paste0(file_path_root, "all_lr_meta.rds"))
get_cellctype_name(pse_data, pse_meta, all_marker_list,paste0(file_path_root, "all_lr.csv"),
                   paste0(file_path_root, "all_lr_old.csv"),paste0(file_path_root, "all_lr_raw.csv"))
>>>>>>> 3e2f6a6 (add large file)
