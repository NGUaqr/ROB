library(igraph)

# 绘制网络，使用度数进行节点颜色填充
plot_network <- function(g, title, file_name) {
  vertex_degree <- igraph::degree(g)
  
  # 生成梯度颜色
  colors <- colorRampPalette(c("yellow", "red"))(max(vertex_degree) + 1)
  vertex_colors <- colors[vertex_degree + 1]  # 根据度数选择颜色
  
  vertex_labels <- paste(V(g)$name, "\n", "degree: ", vertex_degree, sep="")
  
  plot(g, vertex.label = vertex_labels, vertex.size = vertex_degree * 1.5,  # 缩小节点大小
       vertex.color = vertex_colors, edge.arrow.size = 0.3,  # 缩小边箭头大小
       main = title, cex = 0.7)  # 调整标签字体大小
}


# 生成具有梯度变化的网络
generate_gradient_network <- function(num_nodes, avg_degree) {
  g <- sample_pa(n = num_nodes, power = 1, m = avg_degree, directed = FALSE)
  V(g)$name <- paste0("Gene", 1:num_nodes)
  return(g)
}

# 计算网络的拓扑指标
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
# edge_df <- get.data.frame(g, what = "edges")
# write.csv(edge_df, "C:/Users/huihui1126/Desktop/药物预测方法测试/robedges.csv", row.names = FALSE)
# 逐步删除节点并计算网络指标变化 --模拟网络逐渐破坏
analyze_metrics_change_gradual <- function(num_nodes, avg_degree, iterations) {
  par(mfrow=c(5, 6), mar=c(2, 2, 2, 2))  # 调整每个子图的边距，增加显示空间
  metrics_diff_all <- list()
  for (iter in 1:iterations) {
    g.copy <- g <- generate_gradient_network(num_nodes, avg_degree)
    plot_network(g, paste("Original Network"))
    metrics_list <- list()
    # 初始网络指标
    original_metrics <- calculate_metrics(g)
    metrics_list[["Original"]] <- original_metrics
    # 节点按度数排序
    nodes_sorted_by_degree <- order(degree(g), decreasing = TRUE)
    #i = 3
    for (i in seq(1, length(nodes_sorted_by_degree), by = 1)) {
      if (i > length(nodes_sorted_by_degree)) break
      # 获取要删除的节点名称
      nodes_to_remove <- V(g)$name[nodes_sorted_by_degree[i]]  
      # 删除节点
      g.copy <- delete_vertices(g.copy, nodes_to_remove) ### 注意网络信息递归
      plot_network(g.copy, paste("Network after removing", paste(nodes_sorted_by_degree[1:i], collapse = " ")))
      # 计算删除后的网络指标
      current_metrics <- calculate_metrics(g.copy)
      # 计算与原始网络相比的变化
      metrics_diff <- current_metrics - original_metrics
      node_names <- paste(nodes_to_remove, collapse = ", ")
      # 将删除的节点名称作为行名
      metrics_list[[node_names]] <- metrics_diff
    }
    metrics_df <- do.call(rbind, metrics_list)
    rownames(metrics_df) <- seq(0, nrow(metrics_df) - 1)  # 设置行名为删除的节点数目
    # 保存每次迭代的结果
    metrics_diff_all[[iter]] <- metrics_df
  }
  # 计算每个步骤的平均值
  avg_metrics_diff <- Reduce("+", metrics_diff_all) / iterations
  return(avg_metrics_diff)
}

# 同时删除1/多个节点并计算网络指标变化 --模拟药物进攻疾病网络
analyze_metrics_change <- function(num_nodes, avg_degree, iterations,n,multi) {
  par(mfrow=c(5, 6), mar=c(2, 2, 2, 2))  # 调整每个子图的边距，增加显示空间
  metrics_diff_all <- list()
  for (iter in 1:iterations) {
    g.copy <- g <- generate_gradient_network(num_nodes, avg_degree)
    plot_network(g, paste("Original Network"))
    metrics_list <- list()
    # 初始网络指标
    original_metrics <- calculate_metrics(g)
    metrics_list[["Original"]] <- original_metrics
    # 节点按度数排序
    nodes_sorted_by_degree <- order(degree(g), decreasing = TRUE)
    #i = 2
    for (i in seq(1, length(nodes_sorted_by_degree), by = n)) {
      if (i >= length(nodes_sorted_by_degree)) break
      # 获取要删除的节点名称
      nodes_to_remove <- V(g)$name[nodes_sorted_by_degree[i:(i+n-1)]]  
      # 删除节点
      if(multi == T){
        g.copy <- delete_vertices(g, nodes_to_remove)  ### 注意网络更新，每次只在原网络上删除三个
        plot_network(g.copy, paste("Network after removing", paste(nodes_sorted_by_degree[i:(i+n-1)], collapse = " ")))
        }
      if(multi == F){
        g.copy <- delete_vertices(g.copy, nodes_to_remove) ### 注意逐渐删除
        plot_network(g.copy, paste("Network after removing", paste(nodes_sorted_by_degree[1:i], collapse = " ")))
        } 
      # 计算删除后的网络指标
      current_metrics <- calculate_metrics(g.copy)
      # 计算与原始网络相比的变化
      metrics_diff <- current_metrics - original_metrics
      node_names <- paste(nodes_to_remove, collapse = ", ")
      # 将删除的节点名称作为行名
      metrics_list[[node_names]] <- metrics_diff
    }
    metrics_df <- do.call(rbind, metrics_list)
    rownames(metrics_df) <- seq(0, nrow(metrics_df) - 1)  # 设置行名为删除的节点数目
    # 保存每次迭代的结果
    metrics_diff_all[[iter]] <- metrics_df
  }
  # 计算每个步骤的平均值
  avg_metrics_diff <- Reduce("+", metrics_diff_all) / iterations
  return(avg_metrics_diff)
}

dev.off()
# 设置参数
num_nodes <- 30
avg_degree <- 4
iterations <- 10
n <- 3
# 分析多次删除节点后的网络指标变化的均值
avg_metrics_df <- analyze_metrics_change(num_nodes, avg_degree, iterations,n,multi=T) #### 通过改变by和+n来控制每次删除几个节点
# 打印结果并保存到文件
print(avg_metrics_df)
write.csv(avg_metrics_df, "./output/avg_metrics_change.csv", row.names = TRUE)

# 计算网络拓扑指标变化的相关性
cor_matrix <- cor(avg_metrics_df[complete.cases(avg_metrics_df), ] %>% .[-1,])

print(cor_matrix)

# 可视化相关性矩阵
library(gplots)
par(mar=c(12, 9, 3, 3)) #留白：下、左、上、右
#heatmap(cor_matrix, main="Correlation of Topological Metrics Changes", symm=TRUE)
heatmap.2(cor_matrix,
          #main = "Correlation of Topological Metrics Changes",  # 标题
          trace = "none",  # 不显示trace
          density.info = "none",  # 不显示密度信息
          margins = c(10, 10),  # 增加边距
          cexRow = 1.2,  # 调整行标签大小
          cexCol = 1.2,  # 调整列标签大小
          col = bluered(100),  # 设置颜色
          dendrogram = "both",  # 显示行和列的聚类树
          key.title = NA,  # 移除key的标题
          key.xlab = "Correlation",  # key的x轴标签
          keysize = 1.5,  # 调整key的大小
          lmat = rbind(c(4, 3), c(2, 1)),  # 设置布局矩阵，4代表key，3代表key的标签，2代表标题，1代表热图
          lhei = c(1, 5),  # 设置布局高度，第一行是标题和key，第二行是热图
          lwid = c(1.5, 4)  # 设置布局宽度，第一列是dendrogram，第二列是热图和key
)
dev.off()  # 完成绘图并关闭图形设备


### 指标可视化测试
# 第一步：将 rownames 作为列提取出来
metrics_df <- avg_metrics_df[-1, ] 
metrics_df <- data.frame(step = as.numeric(rownames(metrics_df)), metrics_df)

# reshape 成长格式：step, Metric, Change
long_df <- melt(metrics_df, id.vars = "step", variable.name = "Metric", value.name = "Change")

# 第二步：绘图
ggplot(long_df, aes(x = step, y = Change, color = Metric, group = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(title = "网络拓扑指标变化趋势",
       x = "删除节点重要程度-从高到低",
       y = "指标变化值") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "bottom")
ggplot(long_df, aes(x = step, y = Change)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "steelblue", size = 2) +
  facet_wrap(~ Metric, scales = "free_y") +  # 每个指标用一个子图，自由y轴
  labs(title = "网络拓扑指标变化趋势（分面图）",
       x = "删除节点重要程度-从高到低",
       y = "指标变化值") +
  theme_minimal(base_size = 14)

### 重复100次测试
all_runs <- list()
for (run in 1:100) {
  cat("Running", run, "\n")
  avg_metrics_df <- analyze_metrics_change(num_nodes, avg_degree, iterations, n, multi = TRUE)
  
  # 转为 data.frame，确保有列名
  avg_metrics_df <- as.data.frame(avg_metrics_df)
  # ✅ 删除 Original 行
  avg_metrics_df <- avg_metrics_df[rownames(avg_metrics_df) != "0", ]
  # 添加 step 和 run 信息
  avg_metrics_df$Step <- as.numeric(rownames(avg_metrics_df))
  avg_metrics_df$Run <- run
  
  all_runs[[run]] <- avg_metrics_df
}

# 合并所有结果
combined_df <- do.call(rbind, all_runs)

# 转成长格式，正确指定 ID 列
library(reshape2)
long_df <- melt(combined_df,
                id.vars = c("Step", "Run"),
                variable.name = "Metric",
                value.name = "Change")

# 画图
library(ggplot2)
long_df$Step <- as.factor(long_df$Step)
ggplot(long_df, aes(x = Step, y = Change, group = Run)) +
  geom_line(alpha = 0.8, color = "gray") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", size = 1.2, color = "red") +
  facet_wrap(~ Metric, scales = "free_y") +
  labs(title = "100次模拟中网络指标的变化趋势",
       x = "删除节点-每次删除1个-度值由高到低", y = "指标变化值") +
  theme_minimal(base_size = 16)


### 排序方法测试
df <- as.data.frame(avg_metrics_df[-1,])
orderd_df <- df %>% sort_columns(.)
sort_columns(df)
x <- c(5, 4, 3, 2, 1)
rank(x)
