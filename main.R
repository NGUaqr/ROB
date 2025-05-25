# 加载必要的库
library(igraph)
library(httr)
library(jsonlite)
library(biomaRt)
library(tidyverse)
library(STRINGdb)
library(progress)
### 更换成你自己的工作目录
library(this.path)
setwd(print(this.path::this.dir()))
#### 基因源转换函数 ####
gene_convert <- function(genes) {
  failed_genes <- c()
  if (length(genes) > 0) {  
    tryCatch({
        # 使用biomaRt将人源基因转换为鼠源基因
        convert <- getLDS(
          attributes = "hgnc_symbol",    # 要转换的属性，这里是基因符号
          filters = "hgnc_symbol",       # 过滤参数
          mart = human_mart,             # 基因名的种属来源，设置为鼠
          values = genes,          # 要转换的基因集
          attributesL = "mgi_symbol",    # 要同源转换的目标属性，这里还是基因符号
          martL = mouse_mart,            # 要同源转换的目标种属，设置为鼠
          uniqueRows = TRUE
        )
        if (length(convert$MGI.symbol) == 0) {
          failed_genes <<- c(failed_genes, genes)
          message("转换失败的基因列表: ", paste(failed_genes, collapse = ", "))
          return(genes)  # 转换失败，返回原始基因名
        } else {
          return(convert$MGI.symbol)
        }
      }, error = function(e) {
        failed_genes <<- c(failed_genes, genes)
        message("转换失败的基因列表: ", paste(failed_genes, collapse = ", "))
        return(genes)  # 转换失败，返回原始基因名
      })
  }
}


#### 疾病上下调基因集的获取 ####
get_disease_genes <- function(disease_file,
                              N = NULL,
                              LogFC = 0) {
  disease_data <- read.table(disease_file, header = T,row.names = 1)
  if (!is.null(N)) {
    up_genes <- rownames(head(disease_data[order(-disease_data[, 1]), , drop = FALSE], N))
    down_genes <- rownames(head(disease_data[order(disease_data[, 1]), , drop = FALSE], N))}
  if (!is.null(LogFC)) {
    up_genes <- rownames(disease_data[disease_data[,1] > LogFC, ,drop=F])
    down_genes <- rownames(disease_data[disease_data[,1] < LogFC, ,drop=F])}
  up_genes <- gene_convert(up_genes)
  down_genes <- gene_convert(down_genes)
  return(list(up = up_genes, down = down_genes))
}


#### 药物上下调基因集的获取 ####
get_drug_genes <- function(drug_file,
                           N = 500,  # select topn & tail n
                           Corr = 0.95  # select corr >0.95 &<-0.95 
                           ){
  drug_data <- read.table(drug_file, header = FALSE, row.names = 1)
  drug_data <- drug_data[order(drug_data[, 1], decreasing = TRUE), ,drop = F]
  # 判断是否需要选择 topN 或者 tailN
  if (!is.null(N)) {
    up_genes <- rownames(head(drug_data, N))
    down_genes <- rownames(tail(drug_data, N))
  }
  if (!is.null(Corr)) {
    up_genes <- rownames(drug_data[drug_data$V2 > Corr, ,drop=F])
    down_genes <- rownames(drug_data[drug_data$V2 < -Corr, ,drop=F])
   }
  up_genes <- gene_convert(up_genes)
  down_genes <- gene_convert(down_genes)
  return(list(up = up_genes, down = down_genes))
}



#### 网络拓扑指标的计算 ####
calculate_metrics <- function(g) {
  density <- edge_density(g)
  dcentr <- centr_degree(g)$centralization
  clocentr <- centr_clo(g, mode = "all")$centralization
  betwcentr <- centr_betw(g, directed = FALSE)$centralization
  eigencentr <- centr_eigen(g, directed = FALSE)$centralization
  spath <- mean_distance(g, directed = FALSE)
  cc <- transitivity(g, type = "global")
  ad <- mean(degree(g))
  
  metrics <- c(density, dcentr, clocentr, betwcentr, eigencentr, spath, cc, ad)
  names(metrics) <- c("density", "dcentr", "clocentr", "betwcentr", "eigencentr", "spath", "cc", "ad")
  return(metrics)
}



#### 蛋白互作网络的构建 ####
build_ppi_network <- function(ppi_data, 
                              genes, 
                              protein_mapping, 
                              score_threshold = 400,
                              reaction = "direct"  # 选择indirect或direct
                              ) {
  selected_ppi_data <- ppi_data[ppi_data["combined_score"] >= score_threshold, ]
  # 转换基因符号为蛋白名称
  mapped_genes <- subset(protein_mapping, subset = symbol %in% genes) %>% .[["ensp"]] %>% paste0("10090.", .)
  if(reaction == "direct"){
    # 仅保留基因相关的蛋白互作
    selected_ppi_data <- selected_ppi_data[selected_ppi_data$protein1 %in% mapped_genes & selected_ppi_data$protein2 %in% mapped_genes, ]}
  if(reaction == "indirect"){
    # 仅保留基因相关的蛋白互作
    selected_ppi_data <- selected_ppi_data[selected_ppi_data$protein1 %in% mapped_genes | selected_ppi_data$protein2 %in% mapped_genes, ]}
    # 再次将网络矩阵转变成Gene_symbol
  selected_ppi_data$protein1 <- protein_mapping[match(sub("^10090\\.", "", selected_ppi_data$protein1),protein_mapping$ensp),] %>% .[["symbol"]]
  selected_ppi_data$protein2 <- protein_mapping[match(sub("^10090\\.", "", selected_ppi_data$protein2),protein_mapping$ensp),] %>% .[["symbol"]]
  g <- graph_from_data_frame(selected_ppi_data, directed = FALSE)
  return(g)
}



#### 相关性网络的构建 ####
build_corr_network <- function(genes,
                               disease_file,
                               reaction = "indirect"  # 选择indirect或direct
                               ) {
    if (grepl("Macrophages.*early", disease_file)) {
        # 处理 Macrophages+early 文件
        print(paste("Processing", disease_file, "with early network"))
        network_data <- read.table("./network/Macrophages_5edges.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      } else if (grepl("Macrophages.*late", disease_file)) {
        # 处理 Macrophages+late 文件
        print(paste("Processing", disease_file, "with late network"))
        network_data <- read.table("./network/Macrophages_4edges.tsv",header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      }else if (grepl("Fibroblast.*early", disease_file)) {
        # 处理 Fibroblast+early 文件
        print(paste("Processing", disease_file, "with early network"))
        network_data <- read.table("./network/Fibroblast_4edges.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      } else if (grepl("Fibroblast.*late", disease_file)) {
        # 处理 Fibroblast+late 文件
        print(paste("Processing", disease_file, "with late network"))
        network_data <- read.table("./network/Fibroblast_6edges.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      } else if (grepl("Granulocyte.*early", disease_file)) {
        # 处理 Granulocyte+early 文件
        print(paste("Processing", disease_file, "with early network"))
        network_data <- read.table("./network/Granulocyte_2edges.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      } else {
        print(paste("No matching network found for", disease_file))
      }
  if(reaction == "direct"){
    selected_network_data <- network_data[network_data$V1 %in% genes & network_data$V2 %in% genes, ]}
  if(reaction == "indirect"){  
    selected_network_data <- network_data[network_data$V1 %in% genes | network_data$V2 %in% genes, ]}
    g <- graph_from_data_frame(selected_network_data, directed = FALSE)
    return(g)
  }



#### 绘制网络 ####

plot_network <- function(g, title) {
  # 获得degree
  vertex_degree <- igraph::degree(g)
  # 生成梯度颜色
  colors <- colorRampPalette(c("yellow", "red"))(max(vertex_degree) + 1)
  vertex_colors <- colors[vertex_degree + 1]  # 根据度数选择颜色
  
  vertex_labels <- paste(V(g)$name, "\n", "degree: ", vertex_degree, sep="")
  # 绘图
  plot(g, vertex.label = vertex_labels, vertex.size = vertex_degree * 1.5,  # 缩小节点大小
       vertex.color = vertex_colors, edge.arrow.size = 0.3,  # 缩小边箭头大小
       main = title, cex = 0.7)  # 调整标签字体大小
}


####删除节点稳健性矩阵生成 ####
delete_nodes_and_calculate_metrics <- function(network, 
                                               genes_to_remove_list, 
                                               drug_name) {
  metrics_list <- list()
  # 计算删除前的原始网络指标
  original_metrics <- calculate_metrics(network)
  metrics_list[["original"]] <- original_metrics
  # 删除所有指定的节点并计算网络指标
  for (i in seq(1, length(genes_to_remove_list), by = 1)) {
    genes_to_remove <- intersect(genes_to_remove_list[[i]], V(network)$name)
    g_copy <- delete_vertices(network, genes_to_remove)
    plot_network(g_copy, paste("After removing nodes", paste(genes_to_remove, collapse = ", ")))
    metrics <- calculate_metrics(g_copy)
    metrics_list[[sub("\\.txt","",drug_name[i])]] <- metrics
  }
  return(metrics_list)
}


#### 将list转换为数据框 ####
con.list_to_df <- function(metrics_list,diff){
  if(diff == T){
    metrics_df <- do.call(rbind, metrics_list)
    rownames(metrics_df) <- names(metrics_list)
    baseline <- metrics_df[1, ]
    metrics_diff <- sweep(metrics_df, 2, as.numeric(baseline), FUN = "-")
    return(as.data.frame(metrics_diff))
  }
  metrics_df <- do.call(rbind, metrics_list)
  rownames(metrics_df) <- names(metrics_list)
  return(as.data.frame(metrics_df))
}


#### 定义得分的函数 ####
sort_columns <- function(df, 
                         asc_cols = c("ad", "density", "cc"), 
                         desc_cols = c("spath", "eigencentr", "betwcentr")) {
  ## 将格式统一为数据框
  df <- as.data.frame(df)
  ## "clocentr", "dcentr"这两个指标不考虑，若考虑和asc列放在一起
  # 对指定列进行升序排序
  for (col in asc_cols) {
    ## 处理特殊值 
    # 将 NA, NaN, -Inf 替换为最小值
    df[[col]][is.na(df[[col]]) | df[[col]] == -Inf] <- -666
    # 将 Inf 替换为最大值
    df[[col]][df[[col]] == Inf] <- 666
    # 排序并赋值，越小的值排名越top，赋予的分数越大
    ranks <- rank(df[[col]], ties.method = "average")
    df[[col]] <- length(ranks) + 1 - ranks
  }
  # 对指定列进行降序排序
  for (col in desc_cols) {
    # 将 NA, NaN, Inf 替换为最大值
    df[[col]][is.na(df[[col]]) | df[[col]] == Inf] <- 666
    # 将 -Inf 替换为最小值
    df[[col]][df[[col]] == -Inf] <- -666
    # 排序并赋值，越大的值排名越top，赋予的分数越大
    ranks <- rank(df[[col]], ties.method = "average")
    df[[col]] <- ranks
  }
  return(df[c(asc_cols,desc_cols)])
}

#### 定义upup downdown updown downup综合排序函数 ####
combine_metrics <- function(upup_metrics_diff_score, 
                            downdown_metrics_diff_score, 
                            updown_metrics_diff_score, 
                            downup_metrics_diff_score) {
  
  # 初始化综合排序数据框
  final_metrics_df <- upup_metrics_diff_score  # 复制结构
  final_metrics_df[,] <- 0  # 初始化为0
  
  # 对每个药物进行指标的综合计算和选择
  for (drug in rownames(final_metrics_df)) {
    for (col in colnames(final_metrics_df)) {
      # 计算 upup + downdown 和 updown + downup 的值
      positive_effect <- upup_metrics_diff_score[drug, col] + downdown_metrics_diff_score[drug, col]
      negative_effect <- updown_metrics_diff_score[drug, col] + downup_metrics_diff_score[drug, col]
      
      # 选择综合排序
      if (positive_effect > negative_effect) {
        final_metrics_df[drug, col] <- positive_effect
      } else {
        final_metrics_df[drug, col] <- -negative_effect
      }
    }
  }
  
  return(final_metrics_df)
}

#### 定义RRA整合函数 ####
RAA_aggregate_final <- function(final_metrics_df,  # 整合排序的数据框
                                input_file,        # 输入文件名
                                N = NA,            # 默认情况下，N 会被自动计算为所有输入排名列表中的唯一元素的数量
                                method = "RRA") {  # 默认使用 RRA 方法 #其他方法  'mean'
  # 选择聚合方法
  if ( method == "RRA") {
    library(RobustRankAggreg)
    prefix <- "RRA_"
    # 保留原始行名
    original_rownames <- rownames(final_metrics_df)
    # 根据文件名确定处理方式
    if (grepl("positive", input_file, ignore.case = TRUE)) {
      # 对于正效应文件，直接从大到小排名
      final_metrics_df <- as.data.frame(lapply(final_metrics_df, function(col) {
        (rank(-col, ties.method = "min")/length(col))}))
    } else if (grepl("negative", input_file, ignore.case = TRUE)) {
      # 对于负效应文件，加上负号后从大到小排名
      final_metrics_df <- as.data.frame(lapply(final_metrics_df, function(col) {
        (rank(col, ties.method = "min")/length(col))}))
    } else {
      stop("文件名未能识别为 'positive' 或 'negative'")
    }
    
    # 重新赋值行名
    rownames(final_metrics_df) <- original_rownames
    mat <- as.matrix(final_metrics_df)
    # 使用 aggregateRanks 函数进行排名聚合
    aggregated_result <- aggregateRanks(rmat = mat, method = "RRA")
    # 将聚合结果转换为数据框
    aggregated_result$Rank <- rank(aggregated_result$Score, ties.method = "average")
    result_df <- as.data.frame(aggregated_result)
  }  
  if ( method == "mean") {
    prefix <- "Mean_"
    # 计算每行的均值
    row_means <- rowMeans(final_metrics_df, na.rm = TRUE)
    # 根据文件名判断排序方式
    if (grepl("positive", input_file, ignore.case = TRUE)) {
      # 对于正效应文件，均值越大，排名越高（即 rank 越小）
      ranks <- rank(-row_means, ties.method = "min")
    } else if (grepl("negative", input_file, ignore.case = TRUE)) {
      # 对于负效应文件，均值越小，排名越高（即 rank 越小）
      ranks <- rank(row_means, ties.method = "min")
    } else {
      stop("文件名未能识别为 'positive' 或 'negative'")
    }
    # 将排名结果添加为新的列
    final_metrics_df$Rank <- ranks
    result_df <- final_metrics_df["Rank"]
  }
  # 生成输出文件名
  output_file <- paste0(output_dir,prefix, basename(input_file))
  # 输出结果到文件
  write.table(result_df, file = output_file, sep = "\t", col.names = TRUE, quote = FALSE)
  return(result_df)
}

#### 保存结果到文件的函数 ####
save_metrics_to_file <- function(metrics_df,file_prefix, output_dir) {
  file_name <- paste0(output_dir, "/", file_prefix, "_metrics.txt")
  write.table(metrics_df, file = file_name, sep = "\t", col.names = T, quote = FALSE)
}



#### 内嵌数据 ####
# human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# save(human_mart, file = "human_mart.RData")
# mouse_mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
# save(mouse_mart, file = "mouse_mart.RData")
load("human_mart.RData")
load("mouse_mart.RData")
ppi_data <- read.table("./network/10090.protein.links.v12.0.txt", sep=" ",header=TRUE)
protein_mapping <- read.table("./network/Mus_ensp2symbol.txt", sep="\t",quote = "",header=T)



#### 主函数  ####
mainrob <- function(disease_files,    ## 疾病文件
                    drug_files,       ## 药物文件
                    drug_gene_lists = "drug_gene_lists695.RData", ## 最好提供药物靶点list，包含每个药物的上下调基因
                    output_dir,       ## 输出文件
                    ppi_data,         ## ppi网络信息-自带
                    protein_mapping,  ## 蛋白名映射-自带
                    save.diff = T,    ## 是否输出结果
                    PPI = F,          ## 是否构建PPI网络
                    Corr = T,         ## 是否构建相关性网络-需自己提供/单细胞数据构建Biotip
                    Rank = T,         ## 是否排序
                    Combined = T      ## 是否合并upup+downdown 及updown+downup的比较
                    ){   
  ## 进度条设置
  total_steps <- length(disease_files)
  pb <- progress_bar$new(
    format = "Processing [:bar] :percent eta: :eta",
    total = total_steps, clear = FALSE, width = 60
  )
  ## 生成所有药物上下调基因列表
    if(is.null(drug_gene_lists) == T){
      drug_list <- drug_up_list <- drug_down_list <- list()
    for (drug_file in drug_files) {
      drug_name <- tools::file_path_sans_ext(basename(drug_file))  # 获取文件名并去掉扩展名
      drug_list[[drug_name]] <- drug_name
      drug_genes <- get_drug_genes(drug_file,N = 500,Corr = NULL)
      
      # 将提取的基因添加到列表中
      drug_up_list[[drug_name]] <- drug_genes$up
      drug_down_list[[drug_name]] <- drug_genes$down
      # print(drug_up_list)
      # print(drug_down_list)
    }
    # 保存 drug_up_list 和 drug_down_list以便于后续重复调用
    #save(drug_list, drug_up_list, drug_down_list, file = "drug_gene_lists496.RData")
    save(drug_list, drug_up_list, drug_down_list, file = "drug_gene_lists695.RData")}else {
      # 如果 drug_gene_lists 文件存在，直接加载它们
      load("drug_gene_lists695.RData")
      #load("drug_gene_lists.RData")
    }
  #i = 1
  for (i in 1:length(disease_files)) {
    disease_file <- disease_files[i]
    ## 读取疾病上下调基因集
    disease_genes <- get_disease_genes(disease_file)
    ## 构建疾病上下调基因的PPI网络
    if(PPI == TRUE){
      up_network <- build_ppi_network(ppi_data, disease_genes$up, protein_mapping,reaction = "indirect") #indirect
      plot_network(up_network, paste("Orignal"))
      down_network <- build_ppi_network(ppi_data, disease_genes$down, protein_mapping,reaction = "indirect") ##indirect
      plot_network(down_network, paste("Orignal"))
    }
    if(Corr == TRUE){
      up_network <- build_corr_network(disease_genes$up, disease_file, reaction = "indirect")
      down_network <- build_corr_network(disease_genes$down, disease_file, reaction = "indirect")
      # 检查up_network和down_network是否为空
      if (ecount(up_network) > 0 && ecount(down_network) > 0) {
        plot_network(up_network, paste("Original"))
        plot_network(down_network, paste("Original"))
      } else {
        # 如果网络为空边，使用PPI方法构建网络
        if (ecount(up_network) == 0) {
          warning("Correlation network for up-regulated genes has no edges, switching to PPI-based network.")
          up_network <- build_ppi_network(ppi_data, disease_genes$up, protein_mapping, reaction = "indirect")
          plot_network(up_network, paste("PPI-based for Up Genes"))
        }
        if (ecount(down_network) == 0) {
          warning("Correlation network for down-regulated genes has no edges, switching to PPI-based network.")
          down_network <- build_ppi_network(ppi_data, disease_genes$down, protein_mapping, reaction = "indirect")
          plot_network(down_network, paste("PPI-based for Down Genes"))
        }
      }
    }
    ## 删除药物节点并计算网络指标
    # drug_up_list <- drug_up_list$drug_458
    # drug_down_list <- drug_down_list$drug_458
    upup_metrics_df <- delete_nodes_and_calculate_metrics(up_network, drug_up_list, names(drug_list)) %>% con.list_to_df(.,diff=F)
    downdown_metrics_df <- delete_nodes_and_calculate_metrics(down_network, drug_down_list, names(drug_list)) %>% con.list_to_df(.,diff=F)
    updown_metrics_df <- delete_nodes_and_calculate_metrics(up_network, drug_down_list, names(drug_list)) %>% con.list_to_df(.,diff=F)
    downup_metrics_df <- delete_nodes_and_calculate_metrics(down_network, drug_up_list, names(drug_list)) %>% con.list_to_df(.,diff=F)
    ## 删除药物节点并计算网络指标差值
    upup_metrics_diff <- delete_nodes_and_calculate_metrics(up_network, drug_up_list, names(drug_list))  %>% con.list_to_df(.,diff=T)
    downdown_metrics_diff <- delete_nodes_and_calculate_metrics(down_network, drug_down_list, names(drug_list)) %>% con.list_to_df(.,diff=T)
    updown_metrics_diff <- delete_nodes_and_calculate_metrics(up_network, drug_down_list, names(drug_list)) %>% con.list_to_df(.,diff=T)
    downup_metrics_diff <- delete_nodes_and_calculate_metrics(down_network, drug_up_list, names(drug_list)) %>% con.list_to_df(.,diff=T)
    ## 排序-"ad", "density", "cc"升序  "spath", "eigencentr", "betwcentr"降序
    if(Rank == T){   
      upup_metrics_diff_score  <- upup_metrics_diff[-1,] %>% sort_columns(.)
      downdown_metrics_diff_score <- downdown_metrics_diff[-1,]  %>% sort_columns(.)
      updown_metrics_diff_score <- updown_metrics_diff[-1,]  %>% sort_columns(.)
      downup_metrics_diff_score <- downup_metrics_diff[-1,]  %>% sort_columns(.)
    }   
    if(Combined == T){   
      ## 进行综合排序
      final_metrics_df <- combine_metrics(upup_metrics_diff_score, 
                                          downdown_metrics_diff_score, 
                                          updown_metrics_diff_score, 
                                          downup_metrics_diff_score)
      }   
    ## 输出文件
    if(save.diff == T){
      ## 输出网络指标改变矩阵 ##
      output_dir <- "C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/"
      save_metrics_to_file(upup_metrics_diff, paste( basename(disease_file), "upup", sep = "_"), output_dir)
      save_metrics_to_file(downdown_metrics_diff, paste(basename(disease_file), "downdown", sep = "_"), output_dir)
      save_metrics_to_file(updown_metrics_diff, paste( basename(disease_file), "updown", sep = "_"), output_dir)
      save_metrics_to_file(downup_metrics_diff, paste(basename(disease_file), "downup", sep = "_"), output_dir)
    }
    output_dir <- "C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/"
    save_metrics_to_file(final_metrics_df, paste(basename(disease_file), "final", sep = "_"), output_dir)
  ## 更新进度条
    pb$tick()
  } 
  message("Processing completed!")
  #return(upup_metrics_diff)
}







#### 实际DNB操作 ####
disease_files <- list.files("./diseaselist-deg2function/", full.names = TRUE, pattern = "\\.txt$")
drug_files <- list.files("./druglist-pos1/", full.names = TRUE, pattern = "\\.txt$")
output_dir <- "./output/"
disease_file <- disease_files[1]
disease_files <- disease_files[2]
drug_files <- drug_files[1]
disease_files <- disease_files[4:5]
mainrob(disease_files, 
        drug_files,
        drug_gene_lists = "drug_gene_lists695.RData",
        output_dir,
        ppi_data, 
        protein_mapping, 
        save.diff = F,
        PPI =T,
        Corr = F,
        Rank = T,
        Combined = T)

#### RAA或均值整合多指标 ####
## 手动整理好哪些文件是负效应 哪些文件是正效应后 执行 ##
## positive files
input_files <- list.files("C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/positive/", full.names = TRUE, pattern = "\\.txt$")
output_dir <- "C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/"
for (input_file in input_files) {
  # 读取数据
  final_metrics_df <- read.table(input_file, header = TRUE, row.names = 1)
  # 调用函数进行 RAA 处理
  RAA_aggregate_final(final_metrics_df, input_file, output_dir, method = "mean")
}
## negative files
input_files <- list.files("C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/negative/", full.names = TRUE, pattern = "\\.txt$")
output_dir <- "C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/"
for (input_file in input_files) {
  # 读取数据
  final_metrics_df <- read.table(input_file, header = TRUE, row.names = 1)
  # 调用函数进行 RAA 处理
  RAA_aggregate_final(final_metrics_df, input_file, output_dir, method = "mean")
}



#### 整合正负效应函数及可视化 ####
combine_scores <- function(output_dir, 
                           method = "mean", ## "RRA" "mean"
                           annotation_file = NULL, 
                           color_palette = NULL,
                           annotation_colors = NULL,
                           save.png = F) {
  # 获取所有以 .txt 为后缀的文件
  input_files <- list.files(output_dir, full.names = TRUE, pattern = "\\.txt$")
  # 分别处理 RRA 和 Mean 前缀的文件
  rra_files <- input_files[grepl("^RRA_", basename(input_files), ignore.case = TRUE)]
  mean_files <- input_files[grepl("^Mean_", basename(input_files), ignore.case = TRUE)]
  # 合并函数
  combine_two_files <- function(file1, file2, rank_col) {
    df1 <- read.table(file1, header = TRUE, row.names = 1)
    df2 <- read.table(file2, header = TRUE, row.names = 1)
    # 处理行名：去除 .txt 后缀，替换空格和-为下划线，并转换首字母为小写
    clean_rownames <- function(names) {
      names <- gsub("\\.txt$", "", names)  # 去除 .txt 后缀
      names <- gsub("[ -]", "_", names)    # 替换空格和-为下划线
      paste0(toupper(substring(names, 1, 1)), substring(names, 2))  # 首字母转换为大写
      return(names)
    }
    rownames(df1) <- clean_rownames(rownames(df1))
    rownames(df2) <- clean_rownames(rownames(df2))
    # 合并两个数据框，按行名对齐
    combined_df <- merge(df1, df2, by = "row.names", suffixes = c("_1", "_2"))
    rownames(combined_df) <- clean_rownames(combined_df$Row.names)
    combined_df$Row.names <- NULL
    # 生成 Combined_score 列，选择每行两个数据框中较小的 rank
    combined_df$Combined_score <- pmin(combined_df[[paste0(rank_col, "_1")]], 
                                       combined_df[[paste0(rank_col, "_2")]])
    return(combined_df)
  }
  # 处理 RRA 文件
  if (method == "RRA" & length(rra_files) >= 2) {
    for (i in seq(1, length(rra_files), by = 2)) {
      if (i + 1 <= length(rra_files)) {
        file1 <- rra_files[i]
        file2 <- rra_files[i + 1]
        combined_rra_df <- combine_two_files(file1, file2, rank_col = "Rank")
        
        # 生成输出文件名，包括两个文件的名字
        output_file <- paste0("combined_", basename(file1) %>%  sub(".txt_final_metrics.txt", "", .), "_and_", basename(file2)%>%  sub(".txt_final_metrics.txt", "", .))
        # write.table(combined_rra_df, file = output_file, sep = "\t", col.names = TRUE, quote = FALSE)
        # 绘制热图
        draw_heatmap(combined_rra_df, output_file, annotation_file, save.png = F)
      }
    }
  }
  
  # 处理 Mean 文件
  i = 1 
  if (method == "mean" & length(mean_files) >= 2) {
    for (i in seq(1, length(mean_files), by = 2)) {
      if (i + 1 <= length(mean_files)) {
        file1 <- mean_files[i]
        file2 <- mean_files[i + 1]
        data <- combined_mean_df <- combine_two_files(file1, file2, rank_col = "Rank")
        # 生成输出文件名，包括两个文件的名字
        output_file <- paste0("combined_", basename(file1) %>%  sub(".txt_final_metrics.txt", "", .), "_and_", basename(file2)%>%  sub(".txt_final_metrics.txt", "", .))
        # write.table(combined_mean_df, file = output_file, sep = "\t", col.names = TRUE, quote = FALSE)
        # 绘制热图
        draw_heatmap(combined_mean_df, output_file, annotation_file, save.png = F)
      }
    }
  }
}
# 绘制热图函数
# 绘制热图函数
draw_heatmap <- function(data, 
                         output_file, 
                         annotation_file = NULL,
                         save.png =F) {
  # 检查注释文件行名和数据框行名是否一致
  check_annotation <- function(data, annotation_file) {
    annotation <- read.table(annotation_file, header = TRUE, row.names = 1)
    annotation_names <- rownames(annotation)
    data_names <- rownames(data)
    if(!setequal(annotation_names, data_names)){
      message("注释文件的行名与数据框行名不一致，注意检查。")
      non_intersecting <- union(setdiff(annotation_names, data_names), setdiff(data_names, annotation_names))
      message(paste(non_intersecting, collapse = ","))
    }
    annotation_fliter <- annotation[data_names, , drop = FALSE]
    return(annotation_fliter)
  }
  # 加载并检查行注释文件（如果提供）
  if (!is.null(annotation_file)) {
    annotation <- check_annotation(data, annotation_file)
    all(rownames(data) == rownames(annotation)) # 应返回TRUE
  } else {
    annotation <- NULL
  }
  # 绘制热图
  library(ggplot2)
  library(pheatmap)
  # 检查并处理 'Combined_score' 列
  if (!"Combined_score" %in% colnames(data)) {
    colnames(data)[1] <- "Combined_score"
    print("你输入的数据可能未整合，请注意检查")
  }
  data <- data[, "Combined_score", drop = FALSE]
  data_matrix <- data[order(data[, "Combined_score"]), ,drop = FALSE]
  #提取early类药物和其他药物的排名进行比较：
  print(output_file)
  pheatmap(data_matrix,
           #annotation_row = annotation,
           #annotation_colors = ann_colors,
           color=colorRampPalette(c("firebrick3","white","navy"))(100),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = paste("Heatmap of", output_file))
  #提取early类药物和其他药物的排名进行比较：
  early_ranks <- data$Combined_score[annotation$Time == 'early']
  other_ranks <- data$Combined_score[annotation$Time == 'late']
  early_mean_rank <- mean(early_ranks)
  other_mean_rank <- mean(other_ranks)
  relative_percentage <- (early_mean_rank / other_mean_rank) * 100
  print(paste("Early类药物的平均排名得分:", early_mean_rank))
  print(paste("late类药物的平均排名得分:", other_mean_rank))
  print(paste("Early-late相对百分比:", relative_percentage))
  #进行Mann-Whitney U检验：
  library(coin)
  data$Category <- annotation$Time
  w_df <- data[data$Category != 'long_term', ]
  w_df$Category <- as.factor(w_df$Category)
  test_result <- wilcox_test(Combined_score~Category, data = w_df, distribution = "exact")
  p_value <- pvalue(test_result)
  print(paste("Mann-Whitney U检验 p值:", p_value))
  #提取实验证据类药物和无实验证据药物的排名进行比较：
  clinical_ranks <-  data[!is.na(annotation$Pathway),] %>% .[["Combined_score"]]
  NA_ranks <-  data[is.na(annotation$Pathway),] %>% .[["Combined_score"]]
  clinical_mean_ranks <- mean(clinical_ranks)
  NA_mean_ranks <- mean(NA_ranks)
  relative_percentage <- (clinical_mean_ranks / NA_mean_ranks) * 100
  print(paste("实验证据类药物的平均排名得分:", clinical_mean_ranks))
  print(paste("无证据类药物的平均排名得分:", NA_mean_ranks))
  print(paste("实验-无证据相对百分比:", relative_percentage))
  #进行Mann-Whitney U检验：
  data$Category <- annotation$Pathway
  data[!is.na(data$Category),][["Category"]] <- 'clinical'
  data[is.na(data$Category),][["Category"]] <- '无'
  data$Category <- as.factor(data$Category)
  test_result <- wilcox_test(Combined_score~Category, data = data, distribution = "exact")
  p_value <- pvalue(test_result)
  print(paste("Mann-Whitney U检验 p值:", p_value))
  if(save.png == T){
    # 保存热图为PNG文件
    png_filename <- sub("\\.txt$", "_heatmap.png", output_file)
    pheatmap(data_matrix,
             annotation_row = annotation,
             annotation_colors = ann_colors,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             main = paste("Heatmap of", output_file))
  }
}

# 示例调用
output_dir <- "C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/"
annotation_file <- "C:/Users/huihui1126/Desktop/药物预测方法测试/moa.txt"
Pathway_color <- c("MAPK" = "red","PI3K/Akt/Nrf2/HO-1" = "blue","RhoA/ROCK" = "green","Sonic Hedgehog" = "purple",
  "Hippo/YAP" = "orange","NF-KB" = "brown","Wnt/b-catenin" = "pink","TLR4" = "yellow",
  "NLRP3/IL-1b" = "cyan","Sildenafil" = "#5F80B4","TGF-b/SMADs" = "black","NRF2/HO-1" = "violet",
  "PI3K/Akt/mTOR" = "magenta","TGFb1/TAK1,TGF-b2,TGF-b3,RhoA/ROCK" = "navy","NA" = "white")
Time_color <- c("early" = "red", "late" = "blue", "long_term" = "green", "NA" = "white")
ann_colors <- list(Pathway = Pathway_color, Time = Time_color) # 颜色设置
combine_scores(output_dir, 
               method = "mean",
               annotation_file, 
               color_palette = NULL,
               annotation_colors = ann_colors,
               save.png = F)
## 单侧检验
for (i in 1:length(list.files(output_dir,pattern = "\\.txt$"))) {
  input_files <- list.files(output_dir,pattern = "\\.txt$")
  if(T){ ##单边RAA
  rra_files <- input_files[grepl("^RRA_", basename(input_files), ignore.case = TRUE)]
  filename <- paste(output_dir,rra_files[i],sep = '')
  print(filename)
  df<-read.table(file = filename,row.names = 1,sep="\t",header=TRUE)
  df<-df["Rank"]}
  if(F){ ##单边MEAN
    mean_files <- input_files[grepl("^Mean_", basename(input_files), ignore.case = TRUE)]
    filename <- paste(output_dir,mean_files[i],sep = '')
    print(filename)
    df<-read.table(file = filename,row.names = 1,sep="\t",header=TRUE)}

  draw_heatmap(data = df,
               output_file = filename,
               annotation_file = annotation_file,
               save.png = F)
}
dev.off()



### 496数据处理 
### 单边
moa <- read_excel("C:/Users/huihui1126/Desktop/药物预测方法测试/MOA496.xlsx")
moa <- moa[c("中文名","Name")]
i = 1
output_dir <- "C:/Users/huihui1126/Desktop/药物预测方法测试/rob/output/"
for (i in 1:length(list.files(output_dir,pattern = "\\.txt$"))) {
  input_files <- list.files(output_dir,pattern = "\\.txt$")
  annotation_row <- as.data.frame(moa)
  if(F){ ##单边RAA
    rra_files <- input_files[grepl("^RRA_", basename(input_files), ignore.case = TRUE)]
    filename <- paste(output_dir,rra_files[i],sep = '')
    print(filename)
    df<-read.table(file = filename,row.names = 1,sep="\t",header=TRUE)
    df<-df["Rank"]}
  if(T){ ##单边MEAN
    mean_files <- input_files[grepl("^Mean_", basename(input_files), ignore.case = TRUE)]
    filename <- paste(output_dir,mean_files[i],sep = '')
    print(filename)
    df<-read.table(file = filename,row.names = 1,sep="\t",header=TRUE)
    rownames(df) <- gsub("drug_", "", rownames(df))
    # 1. 将行名转换为普通列
    df$compound_id <- rownames(df)
    annotation_row$compound_id <- rownames(annotation_row)
    # 2. 根据 'compound_id' 进行合并
    merged_df <- merge(df, annotation_row, by = "compound_id", all = TRUE)
    merged_df <- na.omit(merged_df)
    rownames(merged_df) <- merged_df$中文名
    merged_df$Score <- rank(merged_df$Rank, ties.method = "average")
    merged_df <- subset(merged_df,Score<=50) 
    merged_df <- merged_df %>% arrange(Score)
    }
  
  pheatmap(as.matrix(merged_df[,"Score", drop = FALSE]),
           scale = "none",
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           show_rownames = TRUE,
           show_colnames = FALSE,
           color=colorRampPalette(c("firebrick3","white","navy"))(10),
           main = basename(filename),
           fontsize = 10)
}
dev.off()














































































































# ################  MCF103和ALL90 相关性检验####
# library(tidyverse)
# library(reshape2)
# library(pheatmap)
# library(purrr)
# # 修改后的函数，将第一列作为行名处理，并处理重复行名
# read_and_merge_folder <- function(folder_path, prefix) {
#   file_names <- list.files(folder_path, pattern = "\\.txt$", full.names = TRUE)
#   file <- file_names[1]
#   # 读取所有文件，并将基因名设为行名
#   data_list <- lapply(file_names, function(file) {
#     df <- read_delim(file, delim = "\t", col_names = F) 
#     df <- df %>%
#       group_by(X1) %>%
#       dplyr::summarise(across(everything(), mean))
#     #rownames(df) <- df$X1
#     df <- df %>% remove_rownames() %>% column_to_rownames("X1")
#     colnames(df) <- paste0(prefix, "_", basename(file))  # 给列名加上前缀
#     return(df)
#   })
#   
#   # 合并所有数据框
#   combined_df <- Reduce(data_list, function(x, y) {
#     merge(x, y, by = "row.names", all = TRUE)
#   })
#   
#   # 确保合并后的行名唯一
#   rownames(combined_df) <- make.unique(combined_df$Row.names)
#   combined_df <- combined_df[ , -1]
#   
#   return(data_list)
# }
# file_names <- list.files("C:/Users/huihui1126/Desktop/cmap_zhangbo/druglist-pos/", pattern = "\\.txt$", full.names = TRUE)
# folder1_path <- "C:/Users/huihui1126/Desktop/cmap_zhangbo/druglist-pos/"
# folder2_path <- "C:/Users/huihui1126/Desktop/cmap_zhangbo/extract"
# 
# data_group1 <- read_and_merge_folder(folder1_path, "90")
# data_group2 <- read_and_merge_folder(folder2_path, "103")
# data_group1_df <- do.call(cbind, data_group1)
# data_group2_df <- do.call(cbind, data_group2)
# 
# 
# merged_data <- merge(data_group1_df, data_group2_df, by = "row.names", all = TRUE)
# merged_data <- merged_data %>%
#   group_by(Row.names) %>%
#   dplyr::summarise(across(everything(), mean))
# #rownames(df) <- df$X1
# merged_data <- merged_data %>% remove_rownames() %>% column_to_rownames("Row.names")
# merged_data <- na.omit(merged_data)
# # 选择列名中包含 "103" 的列
# cols_103 <- grep("103", colnames(merged_data), value = TRUE)
# # 选择列名中包含 "90" 的列
# cols_90 <- grep("90", colnames(merged_data), value = TRUE)
# # 从每个组中按顺序抽取前 10 列
# sampled_cols_103 <- cols_103[1:10]
# sampled_cols_90 <- cols_90[1:10]
# # 合并选择的列
# selected_cols <- c(sampled_cols_103, sampled_cols_90)
# # 创建一个新的数据框，只包含选定的列
# sampled_df <- merged_data[, selected_cols]
# ## 计算相关性矩阵
# correlation_matrix <- cor(sampled_df, use = "pairwise.complete.obs")
# # 绘制相关性热图
# pheatmap(correlation_matrix, cluster_rows = F,cluster_cols = F,
#          clustering_distance_rows = "euclidean",
#          clustering_distance_cols = "euclidean",
#          color = colorRampPalette(c("blue", "white", "red"))(100),
#          main = "Correlation Heatmap")
# 
# 
# 
# 
# 
# ################ 指标相关性检验证 ####################
# library(igraph)
# 
# # 绘制网络，使用度数进行节点颜色填充
# plot_network <- function(g, title, file_name) {
#   vertex_degree <- degree(g)
#   
#   # 生成梯度颜色
#   colors <- colorRampPalette(c("yellow", "red"))(max(vertex_degree) + 1)
#   vertex_colors <- colors[vertex_degree + 1]  # 根据度数选择颜色
#   
#   vertex_labels <- paste(V(g)$name, "\n", "degree: ", vertex_degree, sep="")
#   
#   plot(g, vertex.label = vertex_labels, vertex.size = vertex_degree * 1.5,  # 缩小节点大小
#        vertex.color = vertex_colors, edge.arrow.size = 0.3,  # 缩小边箭头大小
#        main = title, cex = 0.7)  # 调整标签字体大小
# }
# 
# 
# # 生成具有梯度变化的网络
# generate_gradient_network <- function(num_nodes, avg_degree) {
#   g <- sample_pa(n = num_nodes, power = 1, m = avg_degree, directed = FALSE)
#   V(g)$name <- paste0("Gene", 1:num_nodes)
#   return(g)
# }
# 
# # 计算网络的拓扑指标
# calculate_metrics <- function(g) {
#   density <- edge_density(g)
#   dcentr <- centr_degree(g)$centralization
#   clocentr <- centr_clo(g, mode = "all")$centralization
#   betwcentr <- centr_betw(g, directed = FALSE)$centralization
#   eigencentr <- centr_eigen(g, directed = FALSE)$centralization
#   spath <- mean_distance(g, directed = FALSE)
#   cc <- transitivity(g, type = "global")
#   ad <- mean(degree(g))
#   
#   metrics <- c(density, dcentr, clocentr, betwcentr, eigencentr, spath, cc, ad)
#   names(metrics) <- c("density", "dcentr", "clocentr", "betwcentr", "eigencentr", "spath", "cc", "ad")
#   return(metrics)
# }
# 
# # 逐步删除节点并计算网络指标变化 --模拟网络逐渐破坏
# analyze_metrics_change_gradual <- function(num_nodes, avg_degree, iterations) {
#   par(mfrow=c(5, 6), mar=c(2, 2, 2, 2))  # 调整每个子图的边距，增加显示空间
#   metrics_diff_all <- list()
#   for (iter in 1:iterations) {
#     g.copy <- g <- generate_gradient_network(num_nodes, avg_degree)
#     plot_network(g, paste("Original Network"))
#     metrics_list <- list()
#     # 初始网络指标
#     original_metrics <- calculate_metrics(g)
#     metrics_list[["Original"]] <- original_metrics
#     # 节点按度数排序
#     nodes_sorted_by_degree <- order(degree(g), decreasing = TRUE)
#     #i = 3
#     for (i in seq(1, length(nodes_sorted_by_degree), by = 1)) {
#       if (i > length(nodes_sorted_by_degree)) break
#       # 获取要删除的节点名称
#       nodes_to_remove <- V(g)$name[nodes_sorted_by_degree[i]]  
#       # 删除节点
#       g.copy <- delete_vertices(g.copy, nodes_to_remove) ### 注意网络信息递归
#       plot_network(g.copy, paste("Network after removing", paste(nodes_sorted_by_degree[1:i], collapse = " ")))
#       # 计算删除后的网络指标
#       current_metrics <- calculate_metrics(g.copy)
#       # 计算与原始网络相比的变化
#       metrics_diff <- current_metrics - original_metrics
#       node_names <- paste(nodes_to_remove, collapse = ", ")
#       # 将删除的节点名称作为行名
#       metrics_list[[node_names]] <- metrics_diff
#     }
#     metrics_df <- do.call(rbind, metrics_list)
#     rownames(metrics_df) <- seq(0, nrow(metrics_df) - 1)  # 设置行名为删除的节点数目
#     # 保存每次迭代的结果
#     metrics_diff_all[[iter]] <- metrics_df
#   }
#   # 计算每个步骤的平均值
#   avg_metrics_diff <- Reduce("+", metrics_diff_all) / iterations
#   return(avg_metrics_diff)
# }
# 
# # 同时删除1/多个节点并计算网络指标变化 --模拟药物进攻疾病网络
# analyze_metrics_change <- function(num_nodes, avg_degree, iterations,n,multi) {
#   par(mfrow=c(5, 6), mar=c(2, 2, 2, 2))  # 调整每个子图的边距，增加显示空间
#   metrics_diff_all <- list()
#   for (iter in 1:iterations) {
#     g.copy <- g <- generate_gradient_network(num_nodes, avg_degree)
#     plot_network(g, paste("Original Network"))
#     metrics_list <- list()
#     # 初始网络指标
#     original_metrics <- calculate_metrics(g)
#     metrics_list[["Original"]] <- original_metrics
#     # 节点按度数排序
#     nodes_sorted_by_degree <- order(degree(g), decreasing = TRUE)
#     #i = 2
#     for (i in seq(1, length(nodes_sorted_by_degree), by = n)) {
#       if (i >= length(nodes_sorted_by_degree)) break
#       # 获取要删除的节点名称
#       nodes_to_remove <- V(g)$name[nodes_sorted_by_degree[i:(i+n-1)]]  
#       # 删除节点
#       if(multi == T){
#         g.copy <- delete_vertices(g, nodes_to_remove)  ### 注意网络更新，每次只在原网络上删除三个
#         plot_network(g.copy, paste("Network after removing", paste(nodes_sorted_by_degree[i:(i+n-1)], collapse = " ")))
#       }
#       if(multi == F){
#         g.copy <- delete_vertices(g.copy, nodes_to_remove) ### 注意逐渐删除
#         plot_network(g.copy, paste("Network after removing", paste(nodes_sorted_by_degree[1:i], collapse = " ")))
#       } 
#       # 计算删除后的网络指标
#       current_metrics <- calculate_metrics(g.copy)
#       # 计算与原始网络相比的变化
#       metrics_diff <- current_metrics - original_metrics
#       node_names <- paste(nodes_to_remove, collapse = ", ")
#       # 将删除的节点名称作为行名
#       metrics_list[[node_names]] <- metrics_diff
#     }
#     metrics_df <- do.call(rbind, metrics_list)
#     rownames(metrics_df) <- seq(0, nrow(metrics_df) - 1)  # 设置行名为删除的节点数目
#     # 保存每次迭代的结果
#     metrics_diff_all[[iter]] <- metrics_df
#   }
#   # 计算每个步骤的平均值
#   avg_metrics_diff <- Reduce("+", metrics_diff_all) / iterations
#   return(avg_metrics_diff)
# }
# 
# 
# # 设置参数
# num_nodes <- 30
# avg_degree <- 4
# iterations <- 10
# n <- 3
# # 分析多次删除节点后的网络指标变化的均值
# avg_metrics_df <- analyze_metrics_change(num_nodes, avg_degree, iterations,n,multi=T) #### 通过改变by和+n来控制每次删除几个节点
# # 打印结果并保存到文件
# print(avg_metrics_df)
# write.csv(avg_metrics_df, "./output/avg_metrics_change.csv", row.names = TRUE)
# 
# # 计算网络拓扑指标变化的相关性
# cor_matrix <- cor(avg_metrics_df[complete.cases(avg_metrics_df), ] %>% .[-1,])
# 
# print(cor_matrix)
# 
# # 可视化相关性矩阵
# library(gplots)
# par(mar=c(12, 9, 3, 3)) #留白：下、左、上、右
# #heatmap(cor_matrix, main="Correlation of Topological Metrics Changes", symm=TRUE)
# heatmap.2(cor_matrix,
#           #main = "Correlation of Topological Metrics Changes",  # 标题
#           trace = "none",  # 不显示trace
#           density.info = "none",  # 不显示密度信息
#           margins = c(10, 10),  # 增加边距
#           cexRow = 1.2,  # 调整行标签大小
#           cexCol = 1.2,  # 调整列标签大小
#           col = bluered(100),  # 设置颜色
#           dendrogram = "both",  # 显示行和列的聚类树
#           key.title = NA,  # 移除key的标题
#           key.xlab = "Correlation",  # key的x轴标签
#           keysize = 1.5,  # 调整key的大小
#           lmat = rbind(c(4, 3), c(2, 1)),  # 设置布局矩阵，4代表key，3代表key的标签，2代表标题，1代表热图
#           lhei = c(1, 5),  # 设置布局高度，第一行是标题和key，第二行是热图
#           lwid = c(1.5, 4)  # 设置布局宽度，第一列是dendrogram，第二列是热图和key
# )
# dev.off()  # 完成绘图并关闭图形设备
# 
# 
# ### 排序方法测试
# df <- as.data.frame(avg_metrics_df[-1,])
# orderd_df <- df %>% sort_columns(.)
# sort_columns(df)
# 
# 
# 
# 
# 
# #### 如果有没有转换成功的基因按如下操作 ####
# # 1. 找出含有尚未转换（大写）的基因名的列表项
# is_human_list <- sapply(drug_up_list, function(gene_vec) {
#   any(grepl("^[A-Z0-9-]+$", gene_vec[[1]]))  # 修改正则表达式，允许基因名包含 "-"
# })
# false_list_names <- names(drug_up_list)[is_human_list]
# # 2. 对这些包含人源基因的列表项进行转换
# converted_drug_up_list <- drug_up_list  # 先复制一份
# converted_drug_up_list["Carteolol"] 
# converted_drug_up_list["drug_412"]
# converted_drug_up_list["drug_413"]
# converted_drug_up_list[["drug_412"]] <- gene_convert(converted_drug_up_list$drug_412)
# converted_drug_up_list[["drug_413"]] <- gene_convert(converted_drug_up_list$drug_413)
# drug_up_list <- converted_drug_up_list
# # 同理对downlist一样处理
# is_human_list <- sapply(drug_down_list, function(gene_vec) {
#   any(grepl("^[A-Z0-9-]+$", gene_vec[[1]]))  # 修改正则表达式，允许基因名包含 "-"
# })
# false_list_names <- names(drug_down_list)[is_human_list]
# # 2. 对这些包含人源基因的列表项进行转换
# converted_drug_down_list <- drug_down_list # 先复制一份
# # "drug_148" "drug_269" "drug_288" "drug_289" "drug_344" "drug_412"
# converted_drug_down_list[["drug_412"]] 
# converted_drug_down_list[["drug_412"]] <- gene_convert(converted_drug_down_list[["drug_412"]] )
# drug_down_list <- converted_drug_down_list
