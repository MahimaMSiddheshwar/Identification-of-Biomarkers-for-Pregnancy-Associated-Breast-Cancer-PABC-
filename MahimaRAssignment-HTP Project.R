#a. Download ‘GSE31192’ from Gene Expression Omnibus and Importing the data.

  library(GEOquery)
  my_id <- "GSE31192"
  gse <- getGEO(my_id)
  gse <- gse[[1]]  # Select the first dataset from the list
  gse

#b. Print expression levels
  
  pData(gse) ## print the sample information
  fData(gse) ## print the gene annotation
  exprs(gse) ## print the expression data


#c. Log-normalize the data
  
  summary(exprs(gse))


#d. Create box-plot illustrating the log-normalized data
  
  exprs(gse) <- log2(exprs(gse))
  boxplot(exprs(gse),outline=FALSE)


#e. Load the corresponding phenotype data
  
  library(dplyr)
  sampleInfo <- pData(gse)
  sampleInfo
  sampleInfo <- select(sampleInfo, source_name_ch1,characteristics_ch1.1)
  sampleInfo <- rename(sampleInfo,group = source_name_ch1, patient=characteristics_ch1.1)
  sampleInfo


#f. Calculating the correlation between samples and displayed it as a heat map
  
  library(pheatmap)
  corMatrix <- cor(exprs(gse),use="c")
  pheatmap(corMatrix)
  rownames(sampleInfo)


#g. Recreating the heat map from part e including annotations
  
  colnames(corMatrix)
  rownames(sampleInfo) <- colnames(corMatrix)
  pheatmap(corMatrix,
           annotation_col=sampleInfo)
  
  


#h. Create a scatter plot of the data utilizing Principle Component Analysis (PCA) vectors
  
  library(ggplot2)
  library(ggrepel)

  pca <- prcomp(t(exprs(gse)))

  ## Join the PCs to the sample information
  cbind(sampleInfo, pca$x) %>%
    ggplot(aes(x = PC1, y=PC2, col=group,label=paste("Patient", patient))) + geom_point() + geom_text_repel()


#i. Identify a list of Differentially Expressed Genes (DEGs) printing them along with associated p-value and their adjusted p-value

  library(limma)
  design <- model.matrix(~0+sampleInfo$group)
  design
  colnames(design) <- c("Normal","Tumour")
  summary(exprs(gse))

  cutoff <- median(exprs(gse))
  is_expressed <- exprs(gse) > cutoff

  keep <- rowSums(is_expressed) > 2

  table(keep)

  gse <- gse[keep,]
  fit <- lmFit(exprs(gse), design)
  head(fit$coefficients)
  contrasts <- makeContrasts(Tumour - Normal, levels=design)

  fit2 <- contrasts.fit(fit, contrasts)
  fit2 <- eBayes(fit2)
  topTable(fit2)
  topTable(fit2, coef=1)
  
  # Number of rows in the expression matrix represents the total genes
  total_genes <- nrow(exprs(gse))  
  
  # Print the total number of genes
  cat("Total number of genes analyzed:", total_genes, "\n")
  
  # Identifying upregulated and downregulated genes
  results <- topTable(fit2, sort.by = "P", number = Inf)
  upregulated_genes <- results[results$logFC > 0 & results$adj.P.Val < 0.05,]
  downregulated_genes <- results[results$logFC < 0 & results$adj.P.Val < 0.05,]
  
  # Head
  head(upregulated_genes)
  head(downregulated_genes)
  
  # Count
  count_upregulated <- nrow(upregulated_genes)
  count_downregulated <- nrow(downregulated_genes)
  
  # Print
  cat("Number of upregulated genes:", count_upregulated, "\n")
  cat("Number of downregulated genes:", count_downregulated, "\n")
  
  

#j. Create a volcano plot displaying all of the measured genes, highlight the DEGs with a different color

  anno <- fData(gse)
  anno
  anno <- select(anno, `Gene Symbol`, ENTREZ_GENE_ID)
  fit2$genes <- anno
  topTable(fit2)
  full_results <- topTable(fit2, number=Inf)
  full_results <- tibble::rownames_to_column(full_results,"ID")
  p_cutoff <- 0.05
  fc_cutoff <- 1

  full_results %>%
    mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>%
    ggplot(aes(x = logFC, y = B, col=Significant)) + geom_point()


#k. Create an MA plot displaying all of the measured genes, highlight the DEGs with a different color
  
  library(ggrepel)
  p_cutoff <- 0.01
  fc_cutoff <- 1
  topN <- 30

  full_results %>%
    mutate(Significant = adj.P.Val < p_cutoff, abs(logFC) > fc_cutoff ) %>%
    mutate(Rank = 1:n(), Label = ifelse(Rank < topN, Gene.Symbol,"")) %>%
    ggplot(aes(x = logFC, y = B, col=Significant,label=Label)) + geom_point() + geom_text_repel(col="black")
  
  
  

#l. Using the identified DEGs to determine GO and KEGG Pathway Enrichment Analysis

  library(clusterProfiler)
  library(org.Hs.eg.db)
  
  degs <- subset(full_results, adj.P.Val < p_cutoff & abs(logFC) > fc_cutoff)
  
  # GO Enrichment
  ego <- enrichGO(gene          = degs$ENTREZ_GENE_ID,
                  OrgDb         = org.Hs.eg.db,  # Make sure you have this Bioconductor package installed
                  keyType       = "ENTREZID",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.1,
                  qvalueCutoff  = 0.2)
  
  # KEGG Pathway Enrichment
  ekg <- enrichKEGG(gene         = degs$ENTREZ_GENE_ID,
                    organism     = 'hsa',
                    keyType      = "ncbi-geneid",
                    pvalueCutoff = 0.1,
                    qvalueCutoff = 0.2)

  ego_df <- as.data.frame(ego)
  ekg_df <- as.data.frame(ekg)
  
  nrow(degs)
  head(degs)
  
  # View results
  head(ego_df)
  if (nrow(ekg_df) > 0) {
    head(ekg_df)
  } else {
    print("No significant KEGG pathways found.")
  }
  
  # Viewing the first few results of GO enrichment
  head(summary(ego))
  
  # Bar plot for GO terms
  barplot(ego, showCategory=10)  # Show top 10 categories
  
  # Dot plot for GO terms
  dotplot(ego, showCategory=10)
  
  


 
  
  