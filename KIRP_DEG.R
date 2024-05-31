library(easyTCGA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ComplexHeatmap)
library(enrichplot)
data <- read.csv("/Users/shi/Documents/Methods_Code/KIRC_mRNA_label.csv", header=TRUE, row.names=1)
data <- t(data)
expr <- as.matrix(data[-nrow(data), ])
clusters <- data[nrow(data), ]
clusters1 <- data[nrow(data), ]
clusters2 <- data[nrow(data), ]
clusters3 <- data[nrow(data), ]

clusters1[clusters1 == 1] <- 2
group1 <- factor(clusters1)
deg_res1 <- diff_analysis(exprset = expr, 
                          group = group1,
                          is_count = F,
                          save = F
)
deg_limma1 <- deg_res1$deg_limma
selected_rows1 <- deg_limma1[(abs(deg_limma1$logFC) > 1.5) & (deg_limma1$adj.P.Val < 0.05), ]
selected_row_names1 <- rownames(selected_rows1)
selected_row_names1 <- gsub("\\..*", "", selected_row_names1)

clusters2[clusters2 == 2] <- 0
group2 <- factor(clusters2)
deg_res2 <- diff_analysis(exprset = expr, 
                          group = group2,
                          is_count = F,
                          save = F
)
deg_limma2 <- deg_res2$deg_limma
selected_rows2 <- deg_limma2[(abs(deg_limma2$logFC) > 1.5) & (deg_limma2$adj.P.Val < 0.05), ]
selected_row_names2 <- rownames(selected_rows2)
selected_row_names2 <- gsub("\\..*", "", selected_row_names2)

clusters3[clusters3 == 0] <- 1
group3 <- factor(clusters3)
deg_res3 <- diff_analysis(exprset = expr, 
                          group = group3,
                          is_count = F,
                          save = F
)
deg_limma3 <- deg_res3$deg_limma
selected_rows3 <- deg_limma3[(abs(deg_limma3$logFC) > 1.5) & (deg_limma3$adj.P.Val < 0.05), ]
selected_row_names3 <- rownames(selected_rows3)
selected_row_names3 <- gsub("\\..*", "", selected_row_names3)

gene_vector <- selected_row_names1
eg <- bitr(gene_vector, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
genelist <- eg$ENTREZID
genelist <- unique(genelist)
enrich_result <- enrichGO(gene = genelist, OrgDb = org.Hs.eg.db, ont='ALL', keyType = 'ENTREZID')
result <- enrich_result@result

p3 <- cnetplot(edox, circular = TRUE, colorEdge = TRUE,cex_label_category = 0.8,cex_label_gene = 0.7)