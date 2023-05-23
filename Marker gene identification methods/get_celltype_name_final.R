<<<<<<< HEAD
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
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

           


statics_cell_anno <- function(result_file1, result_file2){
total_num_all=0
total_num_1=0
total_num_2=0
total_percent_1=0
total_precent_2 = 0
  result1 = read.csv(result_file1,row.names=1)
  colnames(result1) = c("Class","Subclass","cluster","Marker1","Marker2","Marker3","Marker4","Marker5","cluster_new")
  original_number = nrow(result1)
  result1 = result1[!is.na(result1$Marker1),]
  print(nrow(result1)/original_number)
  result2 = read.csv(result_file2,row.names=1)
  original_number = nrow(result2)
  colnames(result2) = c("Class","Subclass","cluster","Marker1","Marker2","Marker3","Marker4","Marker5","cluster_new")
  result2 = result2[!is.na(result2$Marker1),]
  print(nrow(result2)/original_number)
  data = rbind(result1, result2)
  cell_type_list = list()
  total_number = 0
  data$use_flag = "no"
  for(i in 1:(nrow(data)-1)){
    #print(i)
    #print(total_number)
    
    if(data[i,"use_flag"]=="no"){
      data[i,"use_flag"]="yes"
      current_class = data[i,"Class"]
      current_subclass = data[i,"Subclass"]
      gene_list = c()
      if(!is.na(data[i,"Marker1"])){
        if(data[i,"Marker1"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker1"])
        }
      }
      if(!is.na(data[i,"Marker2"])){
        if(data[i,"Marker2"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker2"])
        }
      }
      if(!is.na(data[i,"Marker3"])){
        if(data[i,"Marker3"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker3"])
        }
      }
      
      while_flag = 0
      while(while_flag ==0){
        while_flag = 1
        for(j in (i+1):nrow(data)){
          if(data[j,"Class"]==current_class && data[j,"Subclass"]==current_subclass){
            if(data[j,"use_flag"]=="no"){
              if(!is.na(data[j,"Marker1"])){
                if(data[j,"Marker1"]!="nan"){
                  if(data[j,"Marker1"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker1"])
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                    while_flag = 0
                  }
                }
              }
              if(!is.na(data[j,"Marker2"])){
                if(data[j,"Marker2"]!="nan"){
                  if(data[j,"Marker2"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker2"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                  }
                }
              }
              if(!is.na(data[j,"Marker3"])){
                if(data[j,"Marker3"]!="nan"){
                  if(data[j,"Marker3"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker3"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                  }
                }
              }
              
              
              
            }
          }
          
        }
      }
      if(length(gene_list)>0){
        total_number = total_number + 1
      }
      
    }
    
  }
  if(data[nrow(data),"use_flag"]=="no"){
    total_number = total_number + 1
  }
  print(total_number)
  print(dim(data))
  total_num_all = total_number
  print((nrow(data)-total_number)/total_number)
  #1
  data = rbind(result1, result1)
  cell_type_list = list()
  total_number = 0
  data$use_flag = "no"
  for(i in 1:(nrow(data)-1)){
    #print(i)
    #print(total_number)
    
    if(data[i,"use_flag"]=="no"){
      data[i,"use_flag"]="yes"
      current_class = data[i,"Class"]
      current_subclass = data[i,"Subclass"]
      gene_list = c()
      if(!is.na(data[i,"Marker1"])){
        if(data[i,"Marker1"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker1"])
        }
      }
      if(!is.na(data[i,"Marker2"])){
        if(data[i,"Marker2"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker2"])
        }
      }
      if(!is.na(data[i,"Marker3"])){
        if(data[i,"Marker3"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker3"])
        }
      }
      
      while_flag = 0
      while(while_flag ==0){
        while_flag = 1
        for(j in (i+1):nrow(data)){
          if(data[j,"Class"]==current_class && data[j,"Subclass"]==current_subclass){
            if(data[j,"use_flag"]=="no"){
              if(!is.na(data[j,"Marker1"])){
                if(data[j,"Marker1"]!="nan"){
                  if(data[j,"Marker1"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker1"])
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                    while_flag = 0
                  }
                }
              }
              if(!is.na(data[j,"Marker2"])){
                if(data[j,"Marker2"]!="nan"){
                  if(data[j,"Marker2"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker2"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                  }
                }
              }
              if(!is.na(data[j,"Marker3"])){
                if(data[j,"Marker3"]!="nan"){
                  if(data[j,"Marker3"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker3"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                  }
                }
              }
              
              
              
            }
          }
          
        }
      }
      if(length(gene_list)>0){
        total_number = total_number + 1
      }
      
    }
    
  }
  if(data[nrow(data),"use_flag"]=="no"){
    total_number = total_number + 1
  }
  print(total_number)
  print(dim(data))
  total_num_1 = total_number
  print((nrow(data)/2-total_number)/total_number)
total_percent_1 = (nrow(data)/2-total_number)/total_number
  #2
  data = rbind(result2, result2)
  cell_type_list = list()
  total_number = 0
  data$use_flag = "no"
  for(i in 1:(nrow(data)-1)){
    #print(i)
    #print(total_number)
    
    if(data[i,"use_flag"]=="no"){
      data[i,"use_flag"]="yes"
      current_class = data[i,"Class"]
      current_subclass = data[i,"Subclass"]
      gene_list = c()
      if(!is.na(data[i,"Marker1"])){
        if(data[i,"Marker1"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker1"])
        }
      }
      if(!is.na(data[i,"Marker2"])){
        if(data[i,"Marker2"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker2"])
        }
      }
      if(!is.na(data[i,"Marker3"])){
        if(data[i,"Marker3"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker3"])
        }
      }
      
      while_flag = 0
      while(while_flag ==0){
        while_flag = 1
        for(j in (i+1):nrow(data)){
          if(data[j,"Class"]==current_class && data[j,"Subclass"]==current_subclass){
            if(data[j,"use_flag"]=="no"){
              if(!is.na(data[j,"Marker1"])){
                if(data[j,"Marker1"]!="nan"){
                  if(data[j,"Marker1"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker1"])
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                    while_flag = 0
                  }
                }
              }
              if(!is.na(data[j,"Marker2"])){
                if(data[j,"Marker2"]!="nan"){
                  if(data[j,"Marker2"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker2"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                  }
                }
              }
              if(!is.na(data[j,"Marker3"])){
                if(data[j,"Marker3"]!="nan"){
                  if(data[j,"Marker3"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker3"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                  }
                }
              }
              
              
              
            }
          }
          
        }
      }
      if(length(gene_list)>0){
        total_number = total_number + 1
      }
      
    }
    
  }
  if(data[nrow(data),"use_flag"]=="no"){
    total_number = total_number + 1
  }
  print(total_number)
  print(dim(data))
  total_num_2 = total_number
  print((nrow(data)/2-total_number)/total_number)
total_percent_2 = (nrow(data)/2-total_number)/total_number
result_final = 0
if((total_num_1+total_num_2)>total_num_all*2){
result_final = (total_num_1+total_num_2 - total_num_all)/total_num_all * (2-(total_num_1+total_num_2 - total_num_all)/total_num_all)
}else{
result_final = (total_num_1+total_num_2 - total_num_all)/total_num_all
}
return(result_final)
=======
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
as_matrix <- function(mat){
  
  tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
  
  row_pos <- mat@i+1
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  val <- mat@x
  
  for (i in seq_along(val)){
    tmp[row_pos[i],col_pos[i]] <- val[i]
  }
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

           


statics_cell_anno <- function(result_file1, result_file2){
total_num_all=0
total_num_1=0
total_num_2=0
total_percent_1=0
total_precent_2 = 0
  result1 = read.csv(result_file1,row.names=1)
  colnames(result1) = c("Class","Subclass","cluster","Marker1","Marker2","Marker3","Marker4","Marker5","cluster_new")
  original_number = nrow(result1)
  result1 = result1[!is.na(result1$Marker1),]
  print(nrow(result1)/original_number)
  result2 = read.csv(result_file2,row.names=1)
  original_number = nrow(result2)
  colnames(result2) = c("Class","Subclass","cluster","Marker1","Marker2","Marker3","Marker4","Marker5","cluster_new")
  result2 = result2[!is.na(result2$Marker1),]
  print(nrow(result2)/original_number)
  data = rbind(result1, result2)
  cell_type_list = list()
  total_number = 0
  data$use_flag = "no"
  for(i in 1:(nrow(data)-1)){
    #print(i)
    #print(total_number)
    
    if(data[i,"use_flag"]=="no"){
      data[i,"use_flag"]="yes"
      current_class = data[i,"Class"]
      current_subclass = data[i,"Subclass"]
      gene_list = c()
      if(!is.na(data[i,"Marker1"])){
        if(data[i,"Marker1"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker1"])
        }
      }
      if(!is.na(data[i,"Marker2"])){
        if(data[i,"Marker2"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker2"])
        }
      }
      if(!is.na(data[i,"Marker3"])){
        if(data[i,"Marker3"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker3"])
        }
      }
      
      while_flag = 0
      while(while_flag ==0){
        while_flag = 1
        for(j in (i+1):nrow(data)){
          if(data[j,"Class"]==current_class && data[j,"Subclass"]==current_subclass){
            if(data[j,"use_flag"]=="no"){
              if(!is.na(data[j,"Marker1"])){
                if(data[j,"Marker1"]!="nan"){
                  if(data[j,"Marker1"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker1"])
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                    while_flag = 0
                  }
                }
              }
              if(!is.na(data[j,"Marker2"])){
                if(data[j,"Marker2"]!="nan"){
                  if(data[j,"Marker2"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker2"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                  }
                }
              }
              if(!is.na(data[j,"Marker3"])){
                if(data[j,"Marker3"]!="nan"){
                  if(data[j,"Marker3"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker3"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                  }
                }
              }
              
              
              
            }
          }
          
        }
      }
      if(length(gene_list)>0){
        total_number = total_number + 1
      }
      
    }
    
  }
  if(data[nrow(data),"use_flag"]=="no"){
    total_number = total_number + 1
  }
  print(total_number)
  print(dim(data))
  total_num_all = total_number
  print((nrow(data)-total_number)/total_number)
  #1
  data = rbind(result1, result1)
  cell_type_list = list()
  total_number = 0
  data$use_flag = "no"
  for(i in 1:(nrow(data)-1)){
    #print(i)
    #print(total_number)
    
    if(data[i,"use_flag"]=="no"){
      data[i,"use_flag"]="yes"
      current_class = data[i,"Class"]
      current_subclass = data[i,"Subclass"]
      gene_list = c()
      if(!is.na(data[i,"Marker1"])){
        if(data[i,"Marker1"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker1"])
        }
      }
      if(!is.na(data[i,"Marker2"])){
        if(data[i,"Marker2"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker2"])
        }
      }
      if(!is.na(data[i,"Marker3"])){
        if(data[i,"Marker3"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker3"])
        }
      }
      
      while_flag = 0
      while(while_flag ==0){
        while_flag = 1
        for(j in (i+1):nrow(data)){
          if(data[j,"Class"]==current_class && data[j,"Subclass"]==current_subclass){
            if(data[j,"use_flag"]=="no"){
              if(!is.na(data[j,"Marker1"])){
                if(data[j,"Marker1"]!="nan"){
                  if(data[j,"Marker1"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker1"])
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                    while_flag = 0
                  }
                }
              }
              if(!is.na(data[j,"Marker2"])){
                if(data[j,"Marker2"]!="nan"){
                  if(data[j,"Marker2"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker2"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                  }
                }
              }
              if(!is.na(data[j,"Marker3"])){
                if(data[j,"Marker3"]!="nan"){
                  if(data[j,"Marker3"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker3"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                  }
                }
              }
              
              
              
            }
          }
          
        }
      }
      if(length(gene_list)>0){
        total_number = total_number + 1
      }
      
    }
    
  }
  if(data[nrow(data),"use_flag"]=="no"){
    total_number = total_number + 1
  }
  print(total_number)
  print(dim(data))
  total_num_1 = total_number
  print((nrow(data)/2-total_number)/total_number)
total_percent_1 = (nrow(data)/2-total_number)/total_number
  #2
  data = rbind(result2, result2)
  cell_type_list = list()
  total_number = 0
  data$use_flag = "no"
  for(i in 1:(nrow(data)-1)){
    #print(i)
    #print(total_number)
    
    if(data[i,"use_flag"]=="no"){
      data[i,"use_flag"]="yes"
      current_class = data[i,"Class"]
      current_subclass = data[i,"Subclass"]
      gene_list = c()
      if(!is.na(data[i,"Marker1"])){
        if(data[i,"Marker1"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker1"])
        }
      }
      if(!is.na(data[i,"Marker2"])){
        if(data[i,"Marker2"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker2"])
        }
      }
      if(!is.na(data[i,"Marker3"])){
        if(data[i,"Marker3"]!="nan"){
          gene_list = c(gene_list, data[i,"Marker3"])
        }
      }
      
      while_flag = 0
      while(while_flag ==0){
        while_flag = 1
        for(j in (i+1):nrow(data)){
          if(data[j,"Class"]==current_class && data[j,"Subclass"]==current_subclass){
            if(data[j,"use_flag"]=="no"){
              if(!is.na(data[j,"Marker1"])){
                if(data[j,"Marker1"]!="nan"){
                  if(data[j,"Marker1"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker1"])
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                    while_flag = 0
                  }
                }
              }
              if(!is.na(data[j,"Marker2"])){
                if(data[j,"Marker2"]!="nan"){
                  if(data[j,"Marker2"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker2"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                    if(!is.na(data[j,"Marker3"])){
                      if(data[j,"Marker3"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker3"])
                      }
                    }
                  }
                }
              }
              if(!is.na(data[j,"Marker3"])){
                if(data[j,"Marker3"]!="nan"){
                  if(data[j,"Marker3"] %in% gene_list){
                    data[j,"use_flag"]="yes"
                    gene_list = c(gene_list, data[j,"Marker3"])
                    while_flag = 0
                    if(!is.na(data[j,"Marker2"])){
                      if(data[j,"Marker2"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker2"])
                      }
                    }
                    if(!is.na(data[j,"Marker1"])){
                      if(data[j,"Marker1"]!="nan"){
                        gene_list = c(gene_list, data[j,"Marker1"])
                      }
                    }
                  }
                }
              }
              
              
              
            }
          }
          
        }
      }
      if(length(gene_list)>0){
        total_number = total_number + 1
      }
      
    }
    
  }
  if(data[nrow(data),"use_flag"]=="no"){
    total_number = total_number + 1
  }
  print(total_number)
  print(dim(data))
  total_num_2 = total_number
  print((nrow(data)/2-total_number)/total_number)
total_percent_2 = (nrow(data)/2-total_number)/total_number
result_final = 0
if((total_num_1+total_num_2)>total_num_all*2){
result_final = (total_num_1+total_num_2 - total_num_all)/total_num_all * (2-(total_num_1+total_num_2 - total_num_all)/total_num_all)
}else{
result_final = (total_num_1+total_num_2 - total_num_all)/total_num_all
}
return(result_final)
>>>>>>> 3e2f6a6 (add large file)
}