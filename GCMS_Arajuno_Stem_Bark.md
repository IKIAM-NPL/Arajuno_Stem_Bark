Untargeted characterization of Amazonian barks extracts plants by Gas
Chromatography - Mass Spectrometry
================
Jefferson Pastuña
2023-08-08

- <a href="#introduction" id="toc-introduction">Introduction</a>
- <a href="#before-to-start" id="toc-before-to-start">Before to start</a>
- <a href="#notame-workflow" id="toc-notame-workflow">Notame workflow</a>
- <a href="#preprocessing" id="toc-preprocessing">Preprocessing</a>
- <a href="#second-pca-and-loading-plot"
  id="toc-second-pca-and-loading-plot">Second PCA and loading plot</a>

## Introduction

This R Script aims to record the procedure given for metabolic profiling
of 4 species of plant used in Amazonia folk medicine. Each step has a
brief explanation, code and graphics.

The workflow used was taken from [“notame”: Workflow for Non-Targeted
LC–MS Metabolic Profiling](https://doi.org/10.3390/metabo10040135).
Which offers a wide variety of functions to perform metabolomic profile
analysis.

## Before to start

The “notame” package accepts as input a feature table that can be
obtained through software such as MZMine, MSDial, among others. In this
case, the feature table was obtained with the help of MZmine. The
(\*.txt) file was fixed to obtain the final feature table input.

## Notame workflow

As a first step for the analysis,“notame” package and other dependency
packages were installed.

``` r
# Notame package installation
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("antonvsdata/notame")

# Dependency packages installation
#install_dependencies
```

Then, a main path and a log system was added to have a record of each
process executed.

``` r
# Notame library call
library(notame)
```

    ## Loading required package: Biobase

    ## Loading required package: BiocGenerics

    ## 
    ## Attaching package: 'BiocGenerics'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     IQR, mad, sd, var, xtabs

    ## The following objects are masked from 'package:base':
    ## 
    ##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
    ##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
    ##     get, grep, grepl, intersect, is.unsorted, lapply, Map, mapply,
    ##     match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
    ##     Position, rank, rbind, Reduce, rownames, sapply, setdiff, sort,
    ##     table, tapply, union, unique, unsplit, which.max, which.min

    ## Welcome to Bioconductor
    ## 
    ##     Vignettes contain introductory material; view with
    ##     'browseVignettes()'. To cite Bioconductor, see
    ##     'citation("Biobase")', and for packages 'citation("pkgname")'.

    ## Loading required package: futile.logger

    ## Loading required package: ggplot2

    ## Loading required package: magrittr

``` r
# Main path
ppath <- "F:/Jefferson-Pastuna/Arajuno_Stem_Bark/"
# Log system
init_log(log_file = paste0(ppath, "Result/GCMS/GCMS_log.txt"))
```

    ## INFO [2023-08-09 15:03:12] Starting logging

Next, the MZmine suitable feature list was imported.

``` r
data <- read_from_excel(file = "Data/GCMS_to_R.xlsx", sheet = 3, 
                        corner_row = 6, corner_column = "G", 
                        split_by = c("Column", "Ion Mode"))
```

Once the data is read, the next step was to create a MetaboSet in order
to obtain a specific R object.

``` r
modes <- construct_metabosets(exprs = data$exprs, 
                              pheno_data = data$pheno_data, 
                              feature_data = data$feature_data,
                              group_col = "Group")
```

We can visualize the raw data in order to inspect the processing
routines.

``` r
# Data extraction
mode_test <- modes$Rxt5_EI
# Boxplot of raw data
raw_bp <- plot_sample_boxplots(mode_test,
                               order_by = "Group",
                               fill_by = "Group")
# PCA of raw data
raw_pca <- plot_pca(mode_test,
                       center = TRUE,
                       shape = "Group",
                       color = "Group")
# Package to plots visualization in a same windows
#if (!requireNamespace("devtools", quietly = TRUE)) {
#  install.packages("devtools")
#}
#devtools::install_github("thomasp85/patchwork")
library(patchwork)
# Plot
raw_pca + raw_bp
```

![](GCMS_Arajuno_Stem_Bark_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

## Preprocessing

The first step of the preprocessing is to change the features with value
equal to 0 to NA.

``` r
# Data extraction
mode <- modes$Rxt5_EI
# Change 0 value to NA
mode <- mark_nas(mode, value = 0)
```

Then, features with low detection rate are first flagged and then could
be removed. The notame package employs two criteria to select this
features. First, is the feature presence in a percentage of QC
injections, and then the feature presence in a percentage within a
sample group or class.

``` r
# Low detection rate
mode <- flag_detection(mode, qc_limit = 0.80, group_limit = 0.75)
```

    ## INFO [2023-08-09 15:03:17] 
    ## 7% of features flagged for low detection rate

``` r
# Some statistics after low detection algorithm
visualizations(mode, prefix = paste0(ppath, "Figure/GCMS/", "Low_Detection"))
```

With these values, features which that were not detected in the 80% of
the QC injections and 75% of sample groups will be flagged as low
detection rate.

The next step for preprocessing correspond to drift correction. The
drift correction can be applied by smoothed cubic spline regression.

``` r
# Drift correction
corrected <- correct_drift(mode)
```

    ## INFO [2023-08-09 15:08:15] 
    ## Starting drift correction at 2023-08-09 15:08:15
    ## INFO [2023-08-09 15:08:38] Drift correction performed at 2023-08-09 15:08:38
    ## INFO [2023-08-09 15:08:45] Inspecting drift correction results 2023-08-09 15:08:45
    ## INFO [2023-08-09 15:09:03] Drift correction results inspected at 2023-08-09 15:09:03
    ## INFO [2023-08-09 15:09:03] 
    ## Drift correction results inspected, report:
    ## Drift_corrected: 93%,  Missing_QCS: 7%

``` r
corrected <- correct_drift(corrected)   # Second correction to improve drift correction
```

    ## INFO [2023-08-09 15:09:03] 
    ## Starting drift correction at 2023-08-09 15:09:03
    ## INFO [2023-08-09 15:09:25] Drift correction performed at 2023-08-09 15:09:25
    ## INFO [2023-08-09 15:09:25] Inspecting drift correction results 2023-08-09 15:09:25
    ## INFO [2023-08-09 15:09:43] Drift correction results inspected at 2023-08-09 15:09:43
    ## INFO [2023-08-09 15:09:43] 
    ## Drift correction results inspected, report:
    ## Drift_corrected: 93%,  Missing_QCS: 7%

``` r
corrected <- correct_drift(corrected)   # Third correction to improve drift correction
```

    ## INFO [2023-08-09 15:09:43] 
    ## Starting drift correction at 2023-08-09 15:09:43
    ## INFO [2023-08-09 15:10:04] Drift correction performed at 2023-08-09 15:10:04
    ## INFO [2023-08-09 15:10:04] Inspecting drift correction results 2023-08-09 15:10:04
    ## INFO [2023-08-09 15:10:22] Drift correction results inspected at 2023-08-09 15:10:22
    ## INFO [2023-08-09 15:10:22] 
    ## Drift correction results inspected, report:
    ## Drift_corrected: 93%,  Missing_QCS: 7%

``` r
# Flag low quality features
corrected <- flag_quality(corrected, condition = "RSD_r < 0.3 & D_ratio_r < 0.5")
```

    ## INFO [2023-08-09 15:10:22] 
    ## 26% of features flagged for low quality

Then we can visualize the data after drift correction.

``` r
# Boxplot
corr_bp <- plot_sample_boxplots(corrected,
                                      order_by = "Group",
                                      fill_by = "Group")
# PCA
corr_pca <- plot_pca(corrected,
                        center = TRUE,
                        shape = "Group",
                        color = "Group") 
# Plot
corr_pca + corr_bp
```

![](GCMS_Arajuno_Stem_Bark_files/figure-gfm/unnamed-chunk-10-1.png)<!-- -->

Contaminant peaks based on the Process Blank will be removed.

``` r
# Removal of contaminants
corrected_no_blank <- flag_contaminants(corrected,
                                        blank_col = "Group",
                                        blank_label = "Blank",
                                        flag_thresh = 0.05,
                                        flag_label = "Contaminant")
# Removal blank group from dataset
corrected_no_blank <- corrected_no_blank[, corrected_no_blank$Group != "Blank"]
# Some statistics after low detection algorithm
visualizations(corrected_no_blank, prefix = paste0(ppath, "Figure/GCMS/", "No_Blank"))
```

The next step removes the QC from the analysis, since they will not be
needed in subsequent treatments.

``` r
corrected_no_qc <- drop_qcs(corrected_no_blank)
```

We can visualize data without QC.

``` r
# Boxplot
no_qc_bp <- plot_sample_boxplots(corrected_no_qc,
                                 order_by = "Group",
                                 fill_by = "Group")
# PCA
no_qc_pca <- plot_pca(corrected_no_qc,
                      center = TRUE,
                      shape = "Group",
                      color = "Group")
# Plot
no_qc_pca + no_qc_bp
```

![](GCMS_Arajuno_Stem_Bark_files/figure-gfm/unnamed-chunk-13-1.png)<!-- -->

The next step is feature clustering. This step helps us reduce the
number of features of the same molecule that were split due to 70 eV EI
ionization behavior.

``` r
clustered <- cluster_features(corrected_no_qc, rt_window = 1/100,
                              all_features = TRUE,
                              corr_thresh = 0.95,
                              d_thresh = 0.8,
                              plotting = TRUE,
                              prefix = paste0(ppath, "Cluster/GCMS/GCMS_"))
```

    ## INFO [2023-08-09 15:10:50] 
    ## Starting feature clustering at 2023-08-09 15:10:50
    ## INFO [2023-08-09 15:10:50] Finding connections between features in Rxt5_EI
    ## [1] 100
    ## [1] 200
    ## [1] 300
    ## [1] 400
    ## [1] 500
    ## [1] 600
    ## [1] 700
    ## [1] 800
    ## [1] 900
    ## [1] 1000
    ## [1] 1100
    ## [1] 1200
    ## [1] 1300
    ## [1] 1400
    ## [1] 1500
    ## [1] 1600
    ## [1] 1700
    ## [1] 1800
    ## [1] 1900
    ## [1] 2000
    ## [1] 2100
    ## [1] 2200
    ## [1] 2300
    ## [1] 2400
    ## [1] 2500
    ## [1] 2600
    ## [1] 2700
    ## [1] 2800
    ## [1] 2900
    ## [1] 3000
    ## [1] 3100
    ## [1] 3200
    ## [1] 3300
    ## [1] 3400
    ## [1] 3500
    ## [1] 3600
    ## [1] 3700
    ## [1] 3800
    ## [1] 3900
    ## [1] 4000
    ## [1] 4100
    ## [1] 4200
    ## [1] 4300
    ## [1] 4400
    ## [1] 4500
    ## [1] 4600
    ## [1] 4700
    ## [1] 4800
    ## [1] 4900
    ## [1] 5000
    ## [1] 5100
    ## [1] 5200
    ## [1] 5300
    ## [1] 5400
    ## INFO [2023-08-09 15:32:12] Found 67065 connections in Rxt5_EI
    ## INFO [2023-08-09 15:32:12] Found 67065 connections
    ## 326 components found
    ## 
    ## Component 100 / 326 
    ## Component 200 / 326 
    ## Component 300 / 326 
    ## 347 components found
    ## 
    ## Component 100 / 347 
    ## Component 200 / 347 
    ## Component 300 / 347 
    ## 171 components found
    ## 
    ## Component 100 / 171 
    ## 75 components found
    ## 
    ## 50 components found
    ## 
    ## 27 components found
    ## 
    ## 3 components found
    ## 
    ## 2 components found
    ## 
    ## 3 components found
    ## 
    ## 4 components found
    ## 
    ## INFO [2023-08-09 15:32:19] Found 607 clusters of 2 or more features, clustering finished at 2023-08-09 15:32:19

    ## [1] "100 / 1008"

    ## [1] "200 / 1008"

    ## [1] "300 / 1008"

    ## [1] "400 / 1008"

    ## [1] "500 / 1008"

    ## [1] "600 / 1008"

    ## [1] "700 / 1008"

    ## [1] "800 / 1008"

    ## [1] "900 / 1008"

    ## [1] "1000 / 1008"

    ## INFO [2023-08-09 15:38:22] Saved cluster plots to: F:/Jefferson-Pastuna/Arajuno_Stem_Bark/Cluster/GCMS/GCMS_

``` r
compressed <- compress_clusters(clustered)
```

    ## INFO [2023-08-09 15:38:22] Clusters compressed, left with 2129 features

We can inspect PCA plot after clustering algorithm.

``` r
# Boxplot
compr_bp <- plot_sample_boxplots(compressed,
                                    order_by = "Group",
                                    fill_by = "Group")
# PCA
compr_pca <- plot_pca(compressed,
                      center = TRUE,
                      shape = "Group",
                      color = "Group")
# Plot
compr_pca + compr_bp
```

![](GCMS_Arajuno_Stem_Bark_files/figure-gfm/unnamed-chunk-15-1.png)<!-- -->

The next step imputes the data.

``` r
# To clean data
set.seed(21)
imputed <- impute_rf(compressed)
```

    ## INFO [2023-08-09 15:38:23] 
    ## Starting random forest imputation at 2023-08-09 15:38:23
    ## INFO [2023-08-09 15:38:24] Out-of-bag error in random forest imputation: 0.195
    ## INFO [2023-08-09 15:38:24] Random forest imputation finished at 2023-08-09 15:38:24

``` r
# To all data
imputed <- impute_rf(imputed, all_features = TRUE)
```

    ## INFO [2023-08-09 15:38:24] 
    ## Starting random forest imputation at 2023-08-09 15:38:24
    ## INFO [2023-08-09 15:39:00] Out-of-bag error in random forest imputation: 0.302
    ## INFO [2023-08-09 15:39:00] Random forest imputation finished at 2023-08-09 15:39:00

We can inspect PCA plot after imputation.

``` r
# Boxplot
imp_bp <- plot_sample_boxplots(imputed,
                               order_by = "Group",
                               fill_by = "Group")
# PCA
imp_pca <- plot_pca(imputed, center = TRUE,
                    shape = "Group",
                    color = "Group")
# Plot
imp_pca + imp_bp
```

![](GCMS_Arajuno_Stem_Bark_files/figure-gfm/unnamed-chunk-17-1.png)<!-- -->

Finally the data is ready to be exported and proceed with the
statistical analysis.

``` r
save(imputed, file = paste0(ppath, "Result/GCMS/Notame_GCMS_output.RData"))
```

# Second PCA and loading plot

Droping flagged features

``` r
# Extract clean data
no_flag <- drop_flagged(imputed)
# Extracting feature height table
peak_height <- exprs(no_flag)
# Extracting Phenotipic data
pheno_data <- no_flag@phenoData@data
```

Preparing data and transposing feature table.

``` r
# Transposing feature height table
transp_table  <- t(peak_height)
# Changing NA to 0 
transp_table[is.na(transp_table)]=0
# Centering and Scaling features
ei_pca <- prcomp(transp_table, center = TRUE, scale. = TRUE)
```

Plotting PCA results.

``` r
# Library to left_join use
library(dplyr)
# PCA scores
scores <- ei_pca$x %>%                   # Get PC coordinates
  data.frame %>%                         # Convert to data frames
  mutate(Sample_ID = rownames(.)) %>%    # Create a new column with the sample names
  left_join(pheno_data )                 # Adding metadata
# PCA plot
ggplot(scores, aes(PC1, PC2, shape = Species, color = Species)) +
  geom_point(size = 3) +
  guides(x=guide_axis(title = "PC1 (35.93 %)"), y=guide_axis(title = "PC2 (25.78 %)")) +
  theme_classic()
```

![](GCMS_Arajuno_Stem_Bark_files/figure-gfm/unnamed-chunk-21-1.png)<!-- -->

``` r
# Save plot
ggsave('Result/GCMS/GCMS_PCA.png', width = 5, height = 4, device='png', dpi="print")
```

Plotting loading results.

``` r
loadings <- ei_pca$rotation %>%           # Extract loadings
  data.frame(Feature_name = rownames(.))  # New column with feat name
```

Creating an artificial table with Feature name and Compound column.

``` r
# Load a metabolite mane table
metab_name <- readxl::read_excel("Data/example_table.xlsx", 2)
# Creating a new small table of the annotated compounds
ei_compouds <- left_join(metab_name, loadings)
# Plotting results
ggplot(loadings, aes(PC1, PC2)) + 
  geom_point(alpha = 0.2) +
  theme_classic() + 
  geom_point(data = ei_compouds, size = 1) +
  ggrepel::geom_label_repel(data = ei_compouds,
                            aes(label = Compound),
                            box.padding = 0.8,
                            label.padding = 0.27,
                            label.r = 0.3,
                            cex = 3) +
  guides(x=guide_axis(title = "PC1 (35.93 %)"), y=guide_axis(title = "PC2 (21.78 %)")) +
  ggsci::scale_color_aaas()
```

![](GCMS_Arajuno_Stem_Bark_files/figure-gfm/unnamed-chunk-23-1.png)<!-- -->

``` r
# Save plot
ggsave('Result/GCMS/GCMS_Loadings.png', width = 5, height = 4, device='png', dpi="print")
```

Finish a record.

``` r
finish_log()
```

    ## INFO [2023-08-09 15:39:02] Finished analysis. Wed Aug  9 15:39:02 2023
    ## Session info:
    ## 
    ## INFO [2023-08-09 15:39:02] R version 4.2.2 (2022-10-31 ucrt)
    ## INFO [2023-08-09 15:39:02] Platform: x86_64-w64-mingw32/x64 (64-bit)
    ## INFO [2023-08-09 15:39:02] Running under: Windows 7 x64 (build 7601) Service Pack 1
    ## INFO [2023-08-09 15:39:02] 
    ## INFO [2023-08-09 15:39:02] Matrix products: default
    ## INFO [2023-08-09 15:39:02] 
    ## INFO [2023-08-09 15:39:02] locale:
    ## INFO [2023-08-09 15:39:02] [1] LC_COLLATE=English_United States.1252 
    ## INFO [2023-08-09 15:39:02] [2] LC_CTYPE=English_United States.1252   
    ## INFO [2023-08-09 15:39:02] [3] LC_MONETARY=English_United States.1252
    ## INFO [2023-08-09 15:39:02] [4] LC_NUMERIC=C                          
    ## INFO [2023-08-09 15:39:02] [5] LC_TIME=English_United States.1252    
    ## INFO [2023-08-09 15:39:02] 
    ## INFO [2023-08-09 15:39:02] attached base packages:
    ## INFO [2023-08-09 15:39:02] [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## INFO [2023-08-09 15:39:02] 
    ## INFO [2023-08-09 15:39:02] other attached packages:
    ## INFO [2023-08-09 15:39:02] [1] dplyr_1.1.0          patchwork_1.1.2.9000 notame_0.2.0        
    ## INFO [2023-08-09 15:39:02] [4] magrittr_2.0.3       ggplot2_3.4.1.9000   futile.logger_1.4.3 
    ## INFO [2023-08-09 15:39:02] [7] Biobase_2.58.0       BiocGenerics_0.44.0 
    ## INFO [2023-08-09 15:39:02] 
    ## INFO [2023-08-09 15:39:02] loaded via a namespace (and not attached):
    ## INFO [2023-08-09 15:39:02]  [1] ggrepel_0.9.2.9999   Rcpp_1.0.10          gert_1.9.2          
    ## INFO [2023-08-09 15:39:02]  [4] tidyr_1.3.0          digest_0.6.31        foreach_1.5.2       
    ## INFO [2023-08-09 15:39:02]  [7] utf8_1.2.3           cellranger_1.1.0     R6_2.5.1            
    ## INFO [2023-08-09 15:39:02] [10] futile.options_1.0.1 sys_3.4.1            evaluate_0.20       
    ## INFO [2023-08-09 15:39:02] [13] highr_0.10           pillar_1.8.1         itertools_0.1-3     
    ## INFO [2023-08-09 15:39:02] [16] rlang_1.0.6          readxl_1.4.2         rstudioapi_0.14     
    ## INFO [2023-08-09 15:39:02] [19] missForest_1.5       rmarkdown_2.20       textshaping_0.3.6   
    ## INFO [2023-08-09 15:39:02] [22] labeling_0.4.2       Rtsne_0.17           igraph_1.4.1.9003   
    ## INFO [2023-08-09 15:39:02] [25] munsell_0.5.0        compiler_4.2.2       xfun_0.37           
    ## INFO [2023-08-09 15:39:02] [28] pkgconfig_2.0.3      askpass_1.1          systemfonts_1.0.4   
    ## INFO [2023-08-09 15:39:02] [31] pcaMethods_1.90.0    htmltools_0.5.4      openssl_2.0.5       
    ## INFO [2023-08-09 15:39:02] [34] tidyselect_1.2.0     tibble_3.1.8         codetools_0.2-18    
    ## INFO [2023-08-09 15:39:02] [37] randomForest_4.7-1.1 fansi_1.0.4          viridisLite_0.4.1   
    ## INFO [2023-08-09 15:39:02] [40] withr_2.5.0          MASS_7.3-58.1        grid_4.2.2          
    ## INFO [2023-08-09 15:39:02] [43] gtable_0.3.1         lifecycle_1.0.3      formatR_1.14        
    ## INFO [2023-08-09 15:39:02] [46] credentials_1.3.2    scales_1.2.1         zip_2.2.2           
    ## INFO [2023-08-09 15:39:02] [49] cli_3.6.0            stringi_1.7.12       farver_2.1.1        
    ## INFO [2023-08-09 15:39:02] [52] fs_1.6.1             doRNG_1.8.6          ggdendro_0.1.23     
    ## INFO [2023-08-09 15:39:02] [55] ragg_1.2.5           generics_0.1.3       vctrs_0.5.2         
    ## INFO [2023-08-09 15:39:02] [58] cowplot_1.1.2        openxlsx_4.2.5.2     ggsci_3.0.0         
    ## INFO [2023-08-09 15:39:02] [61] lambda.r_1.2.4       RColorBrewer_1.1-3   iterators_1.0.14    
    ## INFO [2023-08-09 15:39:02] [64] tools_4.2.2          glue_1.6.2           purrr_1.0.1         
    ## INFO [2023-08-09 15:39:02] [67] rngtools_1.5.2       parallel_4.2.2       fastmap_1.1.0       
    ## INFO [2023-08-09 15:39:02] [70] yaml_2.3.7           colorspace_2.1-0     knitr_1.42          
    ## INFO [2023-08-09 15:39:02] [73] usethis_2.1.6
