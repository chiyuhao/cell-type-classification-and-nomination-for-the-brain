<<<<<<< HEAD
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
library(SingleCellExperiment)
library(SC3)
library(scater)
library(scCCESS)
library(SingleCellExperiment)
library(reticulate)
library(tensorflow)
library(scCCESS)
library(SIMLR)
file_path_root = "/root/datasets/mouseVISP/"
data1 = read.csv(paste0(file_path_root, "data1/data.csv"),row.names=1)
data2 = read.csv(paste0(file_path_root, "data2/data.csv"),row.names=1)



#consensus1
if(file.exists(paste0(file_path_root,"data1/consesus_1_result_1.rds"))){

result_1 = readRDS(paste0(file_path_root,"data1/consesus_1_result_1.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/consesus_1_result_2.rds"))
result_1_data = data.frame(names(result_1$cl),result_1$cl)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2$cl),result_2$cl)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_consensus1.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_consensus1.rds"))
}
#consensus100
if(file.exists(paste0(file_path_root,"data1/consesus_100_result_1.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/consesus_100_result_1.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/consesus_100_result_2.rds"))
result_1_data = data.frame(names(result_1$cl.result$cl),result_1$cl.result$cl)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2$cl.result$cl),result_2$cl.result$cl)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_consensus100.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_consensus100.rds"))
}
#seurat
if(file.exists(paste0(file_path_root,"data1/seurat_result_10.csv"))){
result_1 = read.csv(paste0(file_path_root,"data1/seurat_result_10.csv"),row.names=1)
result_2 = read.csv(paste0(file_path_root,"data2/seurat_result_10.csv"),row.names=1)
result_1_data = data.frame(rownames(result_1),result_1$seurat_clusters)
rownames(result_1_data)=rownames(result_1)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(rownames(result_2),result_2$seurat_clusters)
rownames(result_2_data)=rownames(result_2)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_seurat.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_seurat.rds"))
}


#monocle
if(file.exists(paste0(file_path_root,"data1/monocle3_10.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/monocle3_10.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/monocle3_10.rds"))
result_1_data = data.frame(names(result_1@clusters$UMAP$clusters),result_1@clusters$UMAP$clusters)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2@clusters$UMAP$clusters),result_2@clusters$UMAP$clusters)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_monocle.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_monocle.rds"))
}


#SC3
if(file.exists(paste0(file_path_root,"data1/SC3.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/SC3.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/SC3.rds"))
result_1_data = data.frame(rownames(result_1@colData),result_1@colData$sc3_150_clusters)
rownames(result_1_data) <- result_1_data[,1]
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(rownames(result_2@colData),result_2@colData$sc3_150_clusters)
rownames(result_2_data) <- result_2_data[,1]
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_sc3.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_sc3.rds"))
}
#Kmeans
if(file.exists(paste0(file_path_root,"data1/scCCESS_kmeans.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/scCCESS_kmeans.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/scCCESS_kmeans.rds"))
result_1_data = data.frame(names(result_1),result_1)
rownames(result_1_data) <- result_1_data[,1]
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2),result_2)
rownames(result_2_data) <- result_2_data[,1]
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_kmeans.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_kmeans.rds"))
}

#simlr
if(file.exists(paste0(file_path_root,"data1/scCCESS_simlr.rds"))){
result_1_name = names(readRDS(paste0(file_path_root,"data1/scCCESS_kmeans.rds")))
result_2_name = names(readRDS(paste0(file_path_root,"data2/scCCESS_kmeans.rds")))
result_1 = readRDS(paste0(file_path_root,"data1/scCCESS_simlr.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/scCCESS_simlr.rds"))
result_1_data = data.frame(result_1_name,result_1)
rownames(result_1_data) <- result_1_data[,1]
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(result_2_name,result_2)
rownames(result_2_data) <- result_2_data[,1]
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_simlr.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_simlr.rds"))
=======
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
library(SingleCellExperiment)
library(SC3)
library(scater)
library(scCCESS)
library(SingleCellExperiment)
library(reticulate)
library(tensorflow)
library(scCCESS)
library(SIMLR)
file_path_root = "/root/datasets/mouseVISP/"
data1 = read.csv(paste0(file_path_root, "data1/data.csv"),row.names=1)
data2 = read.csv(paste0(file_path_root, "data2/data.csv"),row.names=1)



#consensus1
if(file.exists(paste0(file_path_root,"data1/consesus_1_result_1.rds"))){

result_1 = readRDS(paste0(file_path_root,"data1/consesus_1_result_1.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/consesus_1_result_2.rds"))
result_1_data = data.frame(names(result_1$cl),result_1$cl)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2$cl),result_2$cl)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_consensus1.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_consensus1.rds"))
}
#consensus100
if(file.exists(paste0(file_path_root,"data1/consesus_100_result_1.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/consesus_100_result_1.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/consesus_100_result_2.rds"))
result_1_data = data.frame(names(result_1$cl.result$cl),result_1$cl.result$cl)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2$cl.result$cl),result_2$cl.result$cl)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_consensus100.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_consensus100.rds"))
}
#seurat
if(file.exists(paste0(file_path_root,"data1/seurat_result_10.csv"))){
result_1 = read.csv(paste0(file_path_root,"data1/seurat_result_10.csv"),row.names=1)
result_2 = read.csv(paste0(file_path_root,"data2/seurat_result_10.csv"),row.names=1)
result_1_data = data.frame(rownames(result_1),result_1$seurat_clusters)
rownames(result_1_data)=rownames(result_1)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(rownames(result_2),result_2$seurat_clusters)
rownames(result_2_data)=rownames(result_2)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_seurat.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_seurat.rds"))
}


#monocle
if(file.exists(paste0(file_path_root,"data1/monocle3_10.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/monocle3_10.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/monocle3_10.rds"))
result_1_data = data.frame(names(result_1@clusters$UMAP$clusters),result_1@clusters$UMAP$clusters)
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2@clusters$UMAP$clusters),result_2@clusters$UMAP$clusters)
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_monocle.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_monocle.rds"))
}


#SC3
if(file.exists(paste0(file_path_root,"data1/SC3.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/SC3.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/SC3.rds"))
result_1_data = data.frame(rownames(result_1@colData),result_1@colData$sc3_150_clusters)
rownames(result_1_data) <- result_1_data[,1]
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(rownames(result_2@colData),result_2@colData$sc3_150_clusters)
rownames(result_2_data) <- result_2_data[,1]
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_sc3.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_sc3.rds"))
}
#Kmeans
if(file.exists(paste0(file_path_root,"data1/scCCESS_kmeans.rds"))){
result_1 = readRDS(paste0(file_path_root,"data1/scCCESS_kmeans.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/scCCESS_kmeans.rds"))
result_1_data = data.frame(names(result_1),result_1)
rownames(result_1_data) <- result_1_data[,1]
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(names(result_2),result_2)
rownames(result_2_data) <- result_2_data[,1]
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_kmeans.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_kmeans.rds"))
}

#simlr
if(file.exists(paste0(file_path_root,"data1/scCCESS_simlr.rds"))){
result_1_name = names(readRDS(paste0(file_path_root,"data1/scCCESS_kmeans.rds")))
result_2_name = names(readRDS(paste0(file_path_root,"data2/scCCESS_kmeans.rds")))
result_1 = readRDS(paste0(file_path_root,"data1/scCCESS_simlr.rds"))
result_2 = readRDS(paste0(file_path_root,"data2/scCCESS_simlr.rds"))
result_1_data = data.frame(result_1_name,result_1)
rownames(result_1_data) <- result_1_data[,1]
colnames(result_1_data)=c("cell","label")
result_2_data = data.frame(result_2_name,result_2)
rownames(result_2_data) <- result_2_data[,1]
colnames(result_2_data)=c("cell","label")
rogue.res1 <- rogue(data1, labels = result_1_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res1, paste0(file_path_root, "data1/rogue_simlr.rds"))
rogue.res2 <- rogue(data2, labels = result_2_data$label,  samples = 1,platform = "UMI", span = 0.6)
saveRDS(rogue.res2, paste0(file_path_root, "data2/rogue_simlr.rds"))
>>>>>>> 3e2f6a6 (add large file)
}