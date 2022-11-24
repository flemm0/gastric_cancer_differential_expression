GI Cancer Gene Expression Data Visualization
================
Flemming Wu

``` r
library(tidyverse)
library(data.table)
library(ggfortify)
library(ggplot2)
library(factoextra)
library(ComplexHeatmap)
library(pvclust)
library(parallel)
library(dendextend)
library(scales)
library(DiagrammeR)
library(gridExtra)
library(tidytext)
library(readxl)
library(circlize)
```

#### Plot Read Counts For all 8 Samples

``` r
samples <- c("SRR585570", "SRR585571", "SRR585572", "SRR585573",
    "SRR585574", "SRR585575", "SRR585576", "SRR585577")

reads <- c(18883105, 16580585, 27917182, 34687136, 32772474,
    34049244, 31929290, 32457417)

groups <- c("normal gastric tissue", "normal gastric tissue",
    "gastric cancer tissue", "gastric cancer tissue", "gastric cancer tissue",
    "gastric cancer cell line", "gastric cancer cell line", "gastric cancer cell line")

df <- data.frame(samples, reads, groups)
# df

df %>%
    ggplot(aes(y = reads, x = samples)) + geom_col(aes(fill = groups)) +
    scale_fill_manual(values = c("#23E7F7", "#2336F7", "#EEF51F")) +
    scale_y_continuous(labels = scales::comma) + geom_hline(yintercept = mean(reads),
    color = "red", linetype = "dotted") + geom_text(aes(label = format(reads,
    big.mark = ",")), position = position_dodge(width = 1), vjust = -0.5) +
    labs(caption = paste(paste("Minimum reads:", format(min(reads),
        big.mark = ","), sep = " "), paste("Average reads:",
        format(mean(reads), big.mark = ","), sep = " "), paste("Maximum reads:",
        format(max(reads), big.mark = ","), sep = " "), sep = "   "),
        x = "Sample", y = "Number of Reads", fill = "Group") +
    theme(plot.caption.position = "plot")
```

<img src="README_files/figure-gfm/pre-alignment metric plot-1.png" style="display: block; margin: auto;" />

#### Draw flowchart describing general bioinformatic workflow

``` r
grViz("digraph flowchart {
      # node definitions with substituted label text
      node [fontname = Arial, shape = rectangle, color = Lavender, style = filled]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
      tab6 [label = '@@6']
      tab7 [label = '@@7']
      tab8 [label = '@@8']
      tab9 [label = '@@9']

      # edge definitions with the node IDs
      tab1 -> tab3;
      tab2 -> tab3;
      tab3 -> tab4; 
      tab3 -> tab5;
      tab4 -> tab6;
      tab4 -> tab7;
      tab6 -> tab8;
      tab7 -> tab9;
      }

      [1]: 'Download read (FASTQ) files from SRA using SRA-Toolkit'
      [2]: 'Download reference genome (FASTA) and annotation (GTF) files from Ensembl'
      [3]: 'Build index with Bowtie2'
      [4]: 'Align reads to whole genome using Tophat2 (BAM)'
      [5]: 'Align reads to whole transcriptome using Tophat2 (BAM)'
      [6]: 'Obtain FPKM values using Cuffdiff with statistics turned off'
      [7]: 'Obtain p-values and log2 fold differences between normal and cancer tissue using Cuffdiff with statistics on'
      [8]: 'Output file: genes.fpkm_tracking'
      [9]: 'Output file: gene_exp.diff'
      ")
```

    ## PhantomJS not found. You can install it with webshot::install_phantomjs(). If it is installed, please make sure the phantomjs executable can be found via the PATH variable.

<div id="htmlwidget-286df22d1a4db59614d9" style="width:768px;height:576px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-286df22d1a4db59614d9">{"x":{"diagram":"digraph flowchart {\n      # node definitions with substituted label text\n      node [fontname = Arial, shape = rectangle, color = Lavender, style = filled]        \n      tab1 [label = \"Download read (FASTQ) files from SRA using SRA-Toolkit\"]\n      tab2 [label = \"Download reference genome (FASTA) and annotation (GTF) files from Ensembl\"]\n      tab3 [label = \"Build index with Bowtie2\"]\n      tab4 [label = \"Align reads to whole genome using Tophat2 (BAM)\"]\n      tab5 [label = \"Align reads to whole transcriptome using Tophat2 (BAM)\"]\n      tab6 [label = \"Obtain FPKM values using Cuffdiff with statistics turned off\"]\n      tab7 [label = \"Obtain p-values and log2 fold differences between normal and cancer tissue using Cuffdiff with statistics on\"]\n      tab8 [label = \"Output file: genes.fpkm_tracking\"]\n      tab9 [label = \"Output file: gene_exp.diff\"]\n\n      # edge definitions with the node IDs\n      tab1 -> tab3;\n      tab2 -> tab3;\n      tab3 -> tab4; \n      tab3 -> tab5;\n      tab4 -> tab6;\n      tab4 -> tab7;\n      tab6 -> tab8;\n      tab7 -> tab9;\n      }","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>

#### Post-Tophat2 Alignment Metrics

``` r
df <- df %>%
    mutate(input = c(37766210, 33161170, 55834364, 69374272,
        65544948, 68098488, 63858580, 64914834), mapped = c(31498668,
        27402515, 45797812, 58072751, 52477113, 55003238, 42375461,
        44831741), pct_mapped = round(mapped/input * 100, 1))

df %>%
    ggplot(aes(y = pct_mapped, x = samples)) + geom_col(aes(fill = groups)) +
    scale_fill_manual(values = c("#23E7F7", "#2336F7", "#EEF51F")) +
    geom_text(aes(label = paste(pct_mapped, "%", sep = " ")),
        position = position_dodge(width = 1), vjust = -0.5) +
    labs(x = "Sample", y = "Percent of Reads Aligned to Genome",
        fill = "Group")
```

<img src="README_files/figure-gfm/alignment metrics-1.png" style="display: block; margin: auto;" />

#### PCA analysis on all genes across all samples

``` r
fpkm <- data.table::fread("./data/genes.fpkm_tracking", sep = "\t")

names(fpkm)[26:29] <- names(fpkm)[26:29] %>% # fix column names that were incorrectly named in Cuffdiff command
  str_replace("ct2", "ct3")


fpkm <- fpkm %>% 
  select(tracking_id, ends_with("_FPKM"))

row.names(fpkm) <- fpkm$tracking_id

fpkm <- fpkm %>% select(ends_with("FPKM"))

fpkm <- t(fpkm)

# normalize everything using log10 conversion
fpkm_norm <- fpkm + 0.1

fpkm_norm <- log10(fpkm_norm)

# create sample column for labeling
fpkm_norm <- as.data.frame(fpkm_norm)
fpkm_norm$samples <- row.names(fpkm_norm)

fpkm_norm.data <- fpkm_norm[-ncol(fpkm_norm)]

fpkm_norm$type <- ifelse(grepl("norm", fpkm_norm$samples), "normal tissue",
                            ifelse(grepl("ct", fpkm_norm$samples), "cancer tissue", "cancer cell line"))
```

``` r
pca <- prcomp(fpkm_norm.data, scale. = FALSE)

pca$x
```

    ##                  PC1       PC2        PC3        PC4         PC5         PC6
    ## norm1_FPKM -65.16371  20.55190 -13.229521   6.388875  -4.3239111   3.9006219
    ## norm2_FPKM -66.66014  14.63853 -14.183265   8.065955  -2.1647517  -2.0809625
    ## ct1_FPKM    36.41180  37.63755  18.874329  -4.732123 -33.4918958 -20.9838262
    ## ct2_FPKM    33.65106  38.53751  20.125099  -1.521727   3.2177481  34.8346724
    ## ct3_FPKM    13.50041  21.53113  17.623186  -5.890797  41.5268135 -18.9785852
    ## ccl1_FPKM   38.13282 -18.68783 -64.021749 -31.834235   0.6476951   1.0735397
    ## ccl2_FPKM  -19.91832 -70.07453  40.320196 -31.700342  -5.2522864   2.8932757
    ## ccl3_FPKM   30.04608 -44.13424  -5.508276  61.224393  -0.1594117  -0.6587357
    ##                     PC7           PC8
    ## norm1_FPKM  29.34359275 -3.305060e-14
    ## norm2_FPKM -29.55169606  4.981894e-13
    ## ct1_FPKM     0.26795485 -1.740182e-13
    ## ct2_FPKM    -3.77078120 -1.042161e-13
    ## ct3_FPKM     2.07326296 -1.115547e-13
    ## ccl1_FPKM    0.07812399 -6.032099e-14
    ## ccl2_FPKM    0.36978174 -2.528532e-13
    ## ccl3_FPKM    1.18976097  2.441909e-13

``` r
autoplot(pca, scale = FALSE, data = fpkm_norm, colour = "type") +
    ggtitle("PCA") + theme_bw()
```

<img src="README_files/figure-gfm/plot pca and edit to scale-1.png" style="display: block; margin: auto;" />

#### Boxplots for all gene FPKM values across all samples

``` r
df <- as.data.frame(fpkm)
df <- df[colSums(df) > 0]  # filter out genes with 0 expression across all samples
df <- df + 0.1
df <- log10(df)
```

``` r
df %>%
    t() %>%
    reshape::melt() %>%
    arrange(X1) %>%
    rename(sample = X2) %>%
    mutate(group = ifelse(grepl("norm", sample), "normal tissue",
        ifelse(grepl("ct", sample), "cancer tissue", "cancer cell line"))) %>%
    ggplot(aes(x = sample, y = value, fill = group)) + geom_boxplot() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
        hjust = 1)) + labs(x = "Sample", y = "log10(FPKM)", fill = "Group") +
    theme_bw()
```

    ## Warning in type.convert.default(X[[i]], ...): 'as.is' should be specified by the
    ## caller; using TRUE

    ## Warning in type.convert.default(X[[i]], ...): 'as.is' should be specified by the
    ## caller; using TRUE

<img src="README_files/figure-gfm/fpkm boxplots of all samples across all genes-1.png" style="display: block; margin: auto;" />

#### Diagram on filtering criteria for Integrated Pathway Analysis

``` r
grViz(diagram = "digraph flowchart {
      # define node aesthetics
      node [fontname = Arial, shape = rect, color = DeepSkyBlue, style = filled, fontcolor = White]        
      tab1 [label = '@@1']
      tab2 [label = '@@2']
      tab3 [label = '@@3']
      tab4 [label = '@@4']
      tab5 [label = '@@5']
# set up node layout
      tab1 -> tab2;
      tab2 -> tab3;
      tab3 -> tab4;
      tab3 -> tab5
      }
      [1]: 'Initial Number of Genes: 62,635'
      [2]: 'Filter p-value < 0.001 and log2 fold change < -5 or > 5'
      [3]: '280 Genes Remaining'
      [4]: '198 Genes Downregulated'
      [5]: '82 Genes Upregulated'
      ")
```

<div id="htmlwidget-16570ea1098781111bd7" style="width:768px;height:576px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-16570ea1098781111bd7">{"x":{"diagram":"digraph flowchart {\n      # define node aesthetics\n      node [fontname = Arial, shape = rect, color = DeepSkyBlue, style = filled, fontcolor = White]        \n      tab1 [label = \"Initial Number of Genes: 62,635\"]\n      tab2 [label = \"Filter p-value < 0.001 and log2 fold change < -5 or > 5\"]\n      tab3 [label = \"280 Genes Remaining\"]\n      tab4 [label = \"198 Genes Downregulated\"]\n      tab5 [label = \"82 Genes Upregulated\"]\n# set up node layout\n      tab1 -> tab2;\n      tab2 -> tab3;\n      tab3 -> tab4;\n      tab3 -> tab5\n      }","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script>

Perform filtering of genes according to the above diagram

``` r
exp <- data.table::fread("./data/gene_exp.diff", sep = "\t")

# the cutoff will be p-val of < 0.001 and log2 fold change greater than 5 or less than -5

exp_filt <- exp %>%
  arrange(`log2(fold_change)`) %>%
  filter(p_value < 0.001) %>%
  filter(`log2(fold_change)` != -Inf, # Inf or -Inf means gene had 0 expression in one of the groups
         `log2(fold_change)` != Inf) %>%
  filter(`log2(fold_change)` <= -5 | `log2(fold_change)` >= 5)
```

Separate up and down regulated genes

``` r
up_genes <- exp_filt %>%
    filter(`log2(fold_change)` > 0) %>%
    select(gene_id, gene)

down_genes <- exp_filt %>%
    filter(`log2(fold_change)` < 0) %>%
    select(gene_id, gene)
```

#### PCA with filtered genes only

``` r
fpkm_filt <- data.table::fread("./data/genes.fpkm_tracking", sep = "\t")

names(fpkm_filt)[26:29] <- names(fpkm_filt)[26:29] %>% # fix column names
  str_replace("ct2", "ct3")


fpkm_filt <- fpkm_filt %>% 
  select(tracking_id, ends_with("_FPKM")) %>%
  filter(tracking_id %in% exp_filt$test_id) %>%
  column_to_rownames("tracking_id")


fpkm_filt_norm <- log10(fpkm_filt + 0.1)          

fpkm_filt_norm <- t(fpkm_filt_norm)

fpkm_filt_norm <- fpkm_filt_norm %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(type = ifelse(grepl("norm", sample), "normal tissue",
                         ifelse(grepl("ct", sample), "cancer tissue", "cancer cell line")))
  

fpkm_filt_norm.data <- fpkm_filt_norm %>% select(-c(sample, type))

pca2 <- prcomp(fpkm_filt_norm.data)

pca2$x
```

    ##             PC1       PC2        PC3        PC4         PC5         PC6
    ## [1,] -18.128702  4.881388  2.1012668 -0.3824233  0.30454770 -0.07984922
    ## [2,] -18.450072  3.743089  1.9017633 -0.4220197  0.06823862  0.26245897
    ## [3,]  12.482702  5.945321 -0.3174978  0.6216570  4.16198547 -1.84216731
    ## [4,]  12.298481  6.106674 -1.1292513 -0.2898049 -0.31707131  3.55491391
    ## [5,]   6.802316  5.170019 -3.6776502  0.2852018 -4.25456629 -1.96712942
    ## [6,]   5.876206 -7.563087  7.1705473  5.9858020 -0.80079993 -0.01144860
    ## [7,]  -6.442397 -9.929611 -8.9600685  1.9811198  1.02344958  0.36192641
    ## [8,]   5.561467 -8.353792  2.9108904 -7.7795327 -0.18578383 -0.27870474
    ##              PC7           PC8
    ## [1,]  1.97923408 -1.355686e-14
    ## [2,] -2.02045228  2.857090e-15
    ## [3,] -0.16016323 -2.346647e-15
    ## [4,]  0.08271395  1.143877e-14
    ## [5,] -0.04064251  2.365295e-15
    ## [6,]  0.03181069 -1.526990e-15
    ## [7,]  0.08614545  1.325329e-15
    ## [8,]  0.04135385 -7.632783e-17

``` r
autoplot(pca2, scale = FALSE, data = fpkm_filt_norm, colour = "type") + ggtitle("PCA on Filtered Genes") + theme_bw()
```

<img src="README_files/figure-gfm/redo PCA with filtered genes only-1.png" style="display: block; margin: auto;" />

Compute Dendrograms

Use pearson correlation for distance and ward.D2 and complete for
clustering method

``` r
rownames(fpkm_filt_norm.data) <- fpkm_filt_norm$sample

# try pearson correlation and euclidean distance for
# distance measures
res.pearson.dist <- get_dist(fpkm_filt_norm.data, method = "pearson")
res.euclidean.dist <- get_dist(fpkm_norm.data, method = "euclidean")

# try ward and complete hierarchical clustering methods
res.eucl.ward.hc <- hclust(d = res.euclidean.dist, method = "ward.D2")
res.eucl.complete.hc <- hclust(d = res.euclidean.dist, method = "complete")
res.pear.ward.hc <- hclust(d = res.pearson.dist, method = "ward.D2")
res.pear.complete.hc <- hclust(d = res.pearson.dist, method = "complete")


p1 <- fviz_dend(res.pear.ward.hc, cex = 0.5, k = 3, color_labels_by_k = TRUE,
    rect = TRUE, rect_fill = TRUE, main = paste("Distance:",
        "pearson", "Method:", "ward.D2", sep = " "))


p2 <- fviz_dend(res.pear.complete.hc, cex = 0.5, k = 3, color_labels_by_k = TRUE,
    rect = TRUE, rect_fill = TRUE, main = paste("Distance:",
        "pearson", "Method:", "complete", sep = " "))

p3 <- fviz_dend(res.eucl.ward.hc, cex = 0.5, k = 3, color_labels_by_k = TRUE,
    rect = TRUE, rect_fill = TRUE, main = paste("Distance:",
        "euclidean", "Method:", "ward.D2", sep = " "))

p4 <- fviz_dend(res.eucl.complete.hc, cex = 0.5, k = 3, color_labels_by_k = TRUE,
    rect = TRUE, rect_fill = TRUE, main = paste("Distance:",
        "euclidean", "Method:", "complete", sep = " "))

grid.arrange(p1, p2, p3, p4)
```

<img src="README_files/figure-gfm/Dendrogram of differentially expressed genes ward method-1.png" style="display: block; margin: auto;" />

It appears that using euclidean distance along with Ward’s method gives
us the clustering that we expect to see given our knowledge of the three
groups.

``` r
fviz_dend(res.eucl.ward.hc, cex = 0.5, k = 3, color_labels_by_k = TRUE,
    rect = TRUE, rect_fill = TRUE, show_labels = TRUE, main = "Hierarchical Clustering of Samples Based on Filtered Genes",
    xlab = "Samples")
```

    ## Warning: `guides(<scale> = FALSE)` is deprecated. Please use `guides(<scale> =
    ## "none")` instead.

<img src="README_files/figure-gfm/visualize p3-1.png" style="display: block; margin: auto;" />

Evaluate Dendrograms

``` r
# create cluster for parallel computing
n.cores <- detectCores() * 3/4
cl <- makePSOCKcluster(n.cores)
clusterSetRNGStream(cl, 1234)

# calculate p-value dendrogram for pearson correlation and
# ward.d2 method
res.pear.ward.pv <- parPvclust(cl, t(fpkm_filt_norm.data), method.hclust = "ward.D2",
    method.dist = "correlation", nboot = 1000)
```

    ## Warning in parPvclust(cl, t(fpkm_filt_norm.data), method.hclust = "ward.D2", : "parPvclust" has been integrated into pvclust (with "parallel" option).
    ## It is available for back compatibility but will be unavailable in the future.

    ## Multiscale bootstrap... Done.

``` r
# calculate p-value dendrogram for pearson correlation and
# complete method
res.pear.complete.pv <- parPvclust(cl, t(fpkm_filt_norm.data),
    method.hclust = "complete", method.dist = "correlation",
    nboot = 1000)
```

    ## Warning in parPvclust(cl, t(fpkm_filt_norm.data), method.hclust = "complete", : "parPvclust" has been integrated into pvclust (with "parallel" option).
    ## It is available for back compatibility but will be unavailable in the future.

    ## Multiscale bootstrap... Done.

``` r
# calculate p-value dendrogram for euclidean dist and
# ward.d2 method
res.eucl.complete.pv <- parPvclust(cl, t(fpkm_filt_norm.data),
    method.hclust = "ward.D2", method.dist = "euclidean", nboot = 1000)
```

    ## Warning in parPvclust(cl, t(fpkm_filt_norm.data), method.hclust = "ward.D2", : "parPvclust" has been integrated into pvclust (with "parallel" option).
    ## It is available for back compatibility but will be unavailable in the future.

    ## Multiscale bootstrap... Done.

``` r
stopCluster(cl)

# plot
p1 <- plot(res.pear.ward.pv, hang = -1, cex = 0.5)
```

<img src="README_files/figure-gfm/p-values for ward clustering-1.png" style="display: block; margin: auto;" />

``` r
p2 <- plot(res.pear.complete.pv, hang = -1, cex = 0.5)
```

<img src="README_files/figure-gfm/p-values for ward clustering-2.png" style="display: block; margin: auto;" />

``` r
p3 <- plot(res.eucl.complete.pv, hang = -1, cex = 0.5)
```

<img src="README_files/figure-gfm/p-values for ward clustering-3.png" style="display: block; margin: auto;" />

``` r
# Clusters with AU >= 95% are considered to be strongly
# supported by data.
```

Above, it appears that using pearson correlation and complete clustering
method clustered ccl3_FPKM poorly. The AU values for using both pearson
and euclidean along with Ward’s method gave stronger AU values.

``` r
dend1 <- as.dendrogram(res.pear.ward.hc)
dend2 <- as.dendrogram(res.eucl.ward.hc)

dendlist(dend1, dend2) %>%
    untangle(method = "step1side") %>%
    tanglegram(main_left = "Pearson Correlation", main_right = "Euclidean Distance",
        main = "Ward.D2", lwd = 1, cex_main = 1, margin_inner = 5)
```

<img src="README_files/figure-gfm/view tanglegram of top clusters-1.png" style="display: block; margin: auto;" />

Using pearson correlation for distance metric instead of euclidean
results in cancer cell line 2 clustering with the normal tissues instead
of with the other cancer cell line samples. I have determined that using
euclidean distance and ward’s method to be the best choices with my
data.

#### Plot Heatmap

``` r
Heatmap(data.matrix(fpkm_filt_norm.data), name = "Expression",
    show_column_names = FALSE, row_split = 3, clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2")
```

<img src="README_files/figure-gfm/heatmap wards-1.png" style="display: block; margin: auto;" />

``` r
boxplot <- rowAnnotation(`log10(fpkm)` = anno_boxplot(data.matrix(fpkm_filt_norm.data)),
    gp = gpar(fill = "#CCCCCC"))

hmap <- Heatmap(data.matrix(fpkm_filt_norm.data), name = "Expression",
    show_column_names = FALSE, row_split = 3, clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2", right_annotation = boxplot)

draw(hmap)
```

<img src="README_files/figure-gfm/annotated heatmap wards-1.png" style="display: block; margin: auto;" />

At this point, I submitted the list of up regulated and down regulated
genes to do a simple integrated pathway analysis (software from QIAGEN).
There were many significant pathways returned, many of the down
regulated ones were involved in immune system function whereas many of
the up regulated pathways have been previously implicated in cancer. The
excel files are located in the IPA-results folder.

Moving forward, I will be looking at individual genes that are involved
in a specific pathway of interest to me, carbohydrate metabolism, which
was down regulated in cancerous samples, compared to normal samples.

#### Heatmap of carbohydrate metabolism pathway.

I will look at the “metabolism of carbohydrates” pathway, which is
downregulated in gi cancer.

``` r
meta <- read_xls("./IPA-results/down-genes_diseases_and_functions.xls",
    col_names = TRUE) %>%
    filter(`Diseases or Functions Annotation` == "Metabolism of carbohydrate") %>%
    unnest_tokens(genes, Molecules) %>%
    mutate(genes = toupper(genes)) %>%
    select(genes)

meta
```

    ## # A tibble: 22 × 1
    ##    genes   
    ##    <chr>   
    ##  1 AKR1B10 
    ##  2 ALDH1A1 
    ##  3 B3GALT5 
    ##  4 B3GAT1  
    ##  5 B3GNT8  
    ##  6 B4GALNT2
    ##  7 CD22    
    ##  8 CEBPA   
    ##  9 CSF1R   
    ## 10 GCNT2   
    ## # … with 12 more rows
    ## # ℹ Use `print(n = ...)` to see more rows

``` r
fpkm <- data.table::fread("./data/genes.fpkm_tracking", sep = "\t")

names(fpkm)[26:29] <- names(fpkm)[26:29] %>% # fix column names
  str_replace("ct2", "ct3")


fpkm <- fpkm %>% 
  select(gene_short_name, ends_with("_FPKM"))

meta <- merge(x = meta,
      y = fpkm, 
      all.x = TRUE,
      all.y = FALSE,
      by.x = "genes",
      by.y = "gene_short_name")


meta <- meta %>% 
  column_to_rownames("genes") %>% 
  t()


meta_norm <- log10(meta + 0.1)
```

``` r
Heatmap(meta_norm, name = "expression", show_column_names = TRUE,
    row_split = 3, clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2")
```

<img src="README_files/figure-gfm/simple heatmap of normalized metabolism pathway fpkms-1.png" style="display: block; margin: auto;" />

``` r
Heatmap(meta_norm, name = "expression", show_column_names = TRUE,
    row_split = 3, clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2") +
    rowAnnotation(`log10(fpkm)` = anno_boxplot(meta_norm)) +
    rowAnnotation(rn = anno_text(rownames(meta_norm)))
```

<img src="README_files/figure-gfm/annotated heatmap-1.png" style="display: block; margin: auto;" />

Using the same clustering method (euclidean distance and ward.D2 method)
it appears that cancer cell line 2 becomes its own cluster, and cancer
cell lines 1 and 3 cluster with the cancer tissue.

``` r
# compare other clustering methods with original to check
# to see if it is still the best method for the selected
# pathway

# original clustering method
ht1 <- Heatmap(meta_norm, name = "ht1", show_column_names = TRUE,
    row_split = 3, clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2")


# pearson correlation with ward.D2 clustering
ht2 <- Heatmap(meta_norm, name = "ht2", show_column_names = TRUE,
    row_split = 3, clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
    clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2",
    column_names_gp = gpar(fontsize = 9))

# euclidean distance with complete clustering

ht3 <- Heatmap(meta_norm, name = "ht3", show_column_names = TRUE,
    row_split = 3, clustering_distance_rows = "euclidean", clustering_distance_columns = "euclidean",
    clustering_method_rows = "complete", clustering_method_columns = "complete",
    column_names_gp = gpar(fontsize = 9))

# pearson correlation with ward.D2 clustering

ht4 <- Heatmap(meta_norm, name = "ht4", show_column_names = TRUE,
    row_split = 3, clustering_distance_rows = "pearson", clustering_distance_columns = "pearson",
    clustering_method_rows = "complete", clustering_method_columns = "complete",
    column_names_gp = gpar(fontsize = 9))

ht1 + ht2 + ht3 + ht4
```

<img src="README_files/figure-gfm/other clustering method testing-1.png" style="display: block; margin: auto;" />

When plotting all 4 heatmaps together with the original heatmap as the
main, the separation is clearest with the original. Individually, ht2
makes cell line 3 its own cluster and clusters cell line 2 and cancer
tissue 3 with normal, ht2 is similar to original heatmap, and ht3 makes
cell line 1 and 3 are their own clusters. I conclude that the original
clustering method is still the best option for this pathway.

#### Linear model test on individual genes for significance

``` r
# add PC1 and PC2 to the data frame for correction

lm_df <- merge(x = pca$x, y = meta_norm, by = "row.names") %>%
    select(-c(4:9)) %>%
    rename(group = Row.names) %>%
    mutate(group = if_else(grepl("norm", group), "normal tissue",
        if_else(grepl("ct", group), "cancer tissue", "cancer cell line")))

lm_df
```

    ##              group       PC1       PC2    AKR1B10     ALDH1A1    B3GALT5
    ## 1 cancer cell line  38.13282 -18.68783 -0.8667345  3.35282094 -1.0000000
    ## 2 cancer cell line -19.91832 -70.07453  2.3123741  1.09610669  0.4440073
    ## 3 cancer cell line  30.04608 -44.13424 -0.6727494 -1.00000000 -1.0000000
    ## 4    cancer tissue  36.41180  37.63755 -0.4087573  0.04874858 -0.7608675
    ## 5    cancer tissue  33.65106  38.53751  0.8380818  0.43979044 -0.7992545
    ## 6    cancer tissue  13.50041  21.53113  0.5618310  1.09281184 -0.7371829
    ## 7    normal tissue -65.16371  20.55190  2.5629337  2.73156056  0.6481111
    ## 8    normal tissue -66.66014  14.63853  2.7410144  2.64429873  0.4187867
    ##       B3GAT1       B3GNT8   B4GALNT2       CD22      CEBPA       CSF1R
    ## 1 -0.1330172 -0.914175892 -0.8510668 -0.8684404  0.4195923 -0.86615847
    ## 2 -1.0000000 -0.102927951 -0.9084027  0.6500655  0.7850621 -0.04919475
    ## 3 -0.8488291 -0.536435468 -0.3672897 -0.8425314  0.1483713 -0.83862840
    ## 4 -0.7116332 -0.804568919 -0.9308675 -0.4277094 -0.6488107 -0.53431384
    ## 5 -0.5791551 -0.567762093 -0.9468305 -0.7295716 -0.8397906 -0.88279999
    ## 6 -0.6272518 -0.674539756 -0.8804159 -0.9545166 -0.5400183 -0.65658944
    ## 7  0.5968466  1.328092018  0.4609159  0.5410423  0.3956653  1.28263344
    ## 8  1.2785638  0.004293853  1.0750978  0.5739282  0.7671596  1.34111388
    ##        GCNT2       GCNT4       GPD1      HNF4G      IL1R2      P2RY6     PIK3CG
    ## 1  1.4636317 -0.06657396 -0.8983979  0.1644540 -0.7225120 -0.9862221 -0.9504510
    ## 2  0.3273875  0.44186945 -0.2332764  1.1656538  1.5976107 -0.9258930 -0.9760585
    ## 3  0.4141759  0.35035464 -0.5252259 -0.4326289 -0.7070213 -0.8695607 -1.0000000
    ## 4 -0.3844508 -0.72266681 -0.8398548 -0.9300050 -0.6265197 -0.9246480 -0.8729915
    ## 5 -0.5577272 -0.39444070 -0.7862884 -0.9007525 -0.1074735 -0.8780891 -0.9273583
    ## 6 -0.2156424 -0.34840801 -0.5251226 -0.1293579  0.1726058 -0.9753528 -0.8799624
    ## 7  1.0234418  1.04904349  0.5816687  0.8488355  1.5503030  0.1006462  0.5091824
    ## 8  1.0750503  1.03196158  0.5389719  0.8673844  1.7064182  0.0431264  0.3212068
    ##        PRKCB      PTPRC     SLC2A4      SMPD3        SST       TFF1
    ## 1 -0.1110756  1.1469214 -0.5204290 -0.8587833 -1.0000000 -1.0000000
    ## 2 -1.0000000 -1.0000000 -0.8653586  0.5994333 -1.0000000  3.1504801
    ## 3 -1.0000000 -0.2571626 -0.7828871 -0.3948998 -1.0000000  0.4843340
    ## 4 -0.9736302 -0.8737423 -0.7393234 -0.7598755 -0.5465272 -1.0000000
    ## 5 -0.8115550 -0.9421784 -0.9698236 -0.9839329  0.3299549 -0.3796318
    ## 6 -0.9725425 -0.9019564 -0.5431362 -0.2696821  0.1842740  1.4602572
    ## 7  0.6943402  1.5685690  0.4202760  1.2415290  1.6823581  3.7051090
    ## 8  0.8581021  1.5758999  0.5540249  1.3146424  1.6113281  3.8476405

``` r
# function to return p-value from lm object
lmp <- function(modelobject) {
    if (class(modelobject) != "lm")
        stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1], f[2], f[3], lower.tail = F)
    attributes(p) <- NULL
    return(p)
}


# lm test on all genes
no_correction <- list()  # initialize empty list to hold p-values for non-corrected lm test
correction <- list()  # empty list to hold p-values for corrected lm test
for (i in names(lm_df[4:length(names(lm_df))])) {
    no_correction[i] <- lmp(lm(get(i) ~ group, lm_df))
    correction[i] <- lmp(lm(get(i) ~ group + PC1 + PC2, lm_df))
}


# store results in data frame
meta_pval <- rbind(as.data.frame(no_correction), as.data.frame(correction)) %>%
    t() %>%
    as.data.frame() %>%
    rename(`Non PC Corrected p-value` = V1, `PC Corrected p-value` = V2) %>%
    format(scientific = FALSE) %>%
    rownames_to_column("gene")

meta_pval <- meta_pval %>%
    mutate(`Non PC Corrected p-value` = as.numeric(`Non PC Corrected p-value`),
        `PC Corrected p-value` = as.numeric(`PC Corrected p-value`))
```

Plot p-values for each gene’s significance in determining difference
between cancer and non cancer

``` r
meta_pval %>%
    reshape::melt() %>%
    # mutate(value = as.numeric(value)) %>%
arrange(gene) %>%
    ggplot(aes(x = forcats::fct_reorder(factor(gene), value),
        y = value)) + geom_bar(aes(fill = variable), stat = "identity",
    position = "dodge", width = 0.5) + scale_fill_manual(values = c("#1B32FA",
    "#FA871B")) + theme(axis.text.x = element_text(angle = 90,
    vjust = 0.5, hjust = 1)) + geom_hline(yintercept = 0.05,
    linetype = "dashed") + labs(y = "p-value", x = "gene")
```

<img src="README_files/figure-gfm/plot p-values-1.png" style="display: block; margin: auto;" />

Plot heatmap of carbohydrate gene expression

``` r
Heatmap(meta_norm, name = "expression", show_column_names = TRUE,
    row_split = 3, clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2") +
    rowAnnotation(`log10(fpkm)` = anno_boxplot(meta_norm)) +
    rowAnnotation(rn = anno_text(rownames(meta_norm)))
```

<img src="README_files/figure-gfm/heatmap of carbohydrate metabolism pathway-1.png" style="display: block; margin: auto;" />

It appears that after correcting for the first two principal components,
all genes except for CD22, B3GNT8, B3GAT1, and ALDH1A1 are significant.

- P2RY6: [promotes tumorigenesis by inhibiting
  apoptosis](https://pubmed.ncbi.nlm.nih.gov/29454075/)

- PIK3CG: [inhibits PI3K-Akt/PKB signaling system responsible for
  tumorigenesis and the progression of colorectal
  cancers](https://pubmed.ncbi.nlm.nih.gov/12473596/)

- **GPD1**: [encodes glycerol-3-phosphate dehydrogenase, which plays a
  crucial role in carbohydrate and lipid metabolism by catalyzing
  conversion of DHAP and NADG to G3P (glycerol-3-phosphate) and
  NAD+.](https://en.wikipedia.org/wiki/Glycerol-3-phosphate_dehydrogenase_1)
  Has been implicated as a tumor suppressor in [breast
  cancer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5731876/) and
  [bladder
  cancer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9284842/)

- CSFR1: [involved in production, differentiation, and function of
  macrophages](https://en.wikipedia.org/wiki/Colony_stimulating_factor_1_receptor)

- **GCNT4**: [predicted to be involved in O-glycan processing and
  carbohydrate metabolic
  process.](https://www.ncbi.nlm.nih.gov/gene/51301). [Found to be
  downregulated in gastric
  cancer](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7457769/)

- PRKCB: [encodes protein involved in B cell activation, apoptosis
  induction, endothelial cell proliferation, and intestinal sugar
  absorption.](https://www.ncbi.nlm.nih.gov/gene/5579) [Previously
  implicated in colorectal
  cancer](https://link.springer.com/article/10.1007/s12253-013-9739-5)

- B4GALNT2: [catalyzes last step in biosynthesis of human Sd(a) antigen
  and Cad
  antigen.](https://www.genecards.org/cgi-bin/carddisp.pl?gene=B4GALNT2)
  [Found to be highly expressed in model triple negative breast
  cancer](https://pubmed.ncbi.nlm.nih.gov/34589428/)

- CEBPA: [encodes transcription factor involved in maturation of certain
  blood cells and believed to act as a tumor
  supressor.](https://medlineplus.gov/genetics/gene/cebpa/) [Implicated
  in acute myeloid
  leukemia](https://www.frontiersin.org/articles/10.3389/fonc.2022.806137/full)

- SST: [regulator of endocrine system through interactions with
  pituitary growth hormone, thyroid stimulating hormone, and most
  hormones of GI tract. Also affects rates of neurotransmission in the
  CNS and proliferation of normal and tumorigenic
  cells](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SST)

- **SLC2A4**: [encodes member of K+ dependent Na/Ca exchanger protein
  family](https://www.genecards.org/cgi-bin/carddisp.pl?gene=SLC24A4)
  Not many articles on its implications in cancer.

- PTPRC: [known to be signaling molecules that regulate a variety of
  cellular processes including cell growth, differentiation, mitosis,
  and oncogenic transformation](https://www.ncbi.nlm.nih.gov/gene/5788)

- GCNT2: [encodes enzyme responsible for formation of blood group I
  antigen.](https://en.wikipedia.org/wiki/GCNT2) [Hypomethylated in
  colorectal cancer](https://pubmed.ncbi.nlm.nih.gov/25750292/)

I have selected three genes of interest to me: GPD1, GCNT4, and SLC24A
for further analysis.

``` r
summary(lm(GPD1 ~ group, lm_df))
```

    ## 
    ## Call:
    ## lm(formula = GPD1 ~ group, data = lm_df)
    ## 
    ## Residuals:
    ##        1        2        3        4        5        6        7        8 
    ## -0.34610  0.31902  0.02707 -0.12277 -0.06920  0.19197  0.02135 -0.02135 
    ## 
    ## Coefficients:
    ##                    Estimate Std. Error t value Pr(>|t|)   
    ## (Intercept)         -0.5523     0.1366  -4.043  0.00989 **
    ## groupcancer tissue  -0.1648     0.1932  -0.853  0.43262   
    ## groupnormal tissue   1.1126     0.2160   5.151  0.00361 **
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.2366 on 5 degrees of freedom
    ## Multiple R-squared:  0.8863, Adjusted R-squared:  0.8409 
    ## F-statistic:  19.5 on 2 and 5 DF,  p-value: 0.004355

``` r
p1 <- lm_df %>%
    select(group, GPD1) %>%
    group_by(group) %>%
    summarise(mean = mean(GPD1), sd = sd(GPD1)) %>%
    ggplot(aes(x = forcats::fct_reorder(factor(group), mean),
        y = mean, fill = group)) + geom_col() + geom_errorbar(aes(x = group,
    ymin = mean - sd, ymax = mean + sd), width = 0.4, colour = "orange",
    alpha = 0.9, size = 1.3) + geom_text(aes(label = paste("Average log10(fpkm):",
    round(mean, 2), sep = " "), hjust = 0.5, vjust = ifelse(mean <
    0, 5.5, -2.5)), size = 3) + geom_hline(yintercept = 0) +
    ylim(-2, 2) + labs(x = "", y = "Average log(10) FPKM", title = "GPD1") +
    annotate("text", label = paste("p-value:", round(lmp(lm(GPD1 ~
        group, lm_df)), 6), sep = " "), x = 2, y = 2)

p2 <- lm_df %>%
    select(group, SLC2A4) %>%
    group_by(group) %>%
    summarise(mean = mean(SLC2A4), sd = sd(SLC2A4)) %>%
    ggplot(aes(x = forcats::fct_reorder(factor(group), mean),
        y = mean, fill = group)) + geom_col() + geom_errorbar(aes(x = group,
    ymin = mean - sd, ymax = mean + sd), width = 0.4, colour = "orange",
    alpha = 0.9, size = 1.3) + geom_text(aes(label = paste("Average log10(fpkm):",
    round(mean, 2), sep = " "), hjust = 0.5, vjust = ifelse(mean <
    0, 5.5, -2.5)), size = 3) + geom_hline(yintercept = 0) +
    ylim(-2, 2) + labs(x = "", y = "Average log(10) FPKM", title = "SLC2A4") +
    annotate("text", label = paste("p-value:", round(lmp(lm(SLC2A4 ~
        group, lm_df)), 6), sep = " "), x = 2, y = 2)


grid.arrange(p1, p2, ncol = 1)
```

<img src="README_files/figure-gfm/barplots of individual genes-1.png" style="display: block; margin: auto;" />
