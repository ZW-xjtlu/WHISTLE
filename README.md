An Introduction to *WHISTLE*
================
2020-06-12

Package Installation
====================

To install WHISTLE from Github, use the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges","GenomicFeatures",
                       "BSgenome","IRanges",
                       "S4Vectors","GenomicScores",
                       "AnnotationDbi","Biostrings",
                       "TxDb.Hsapiens.UCSC.hg19.knownGene"))

if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("ZW-xjtlu/WHISTLE", build_vignettes = TRUE)
```

If you want to annotate sequence and evolutionary conservation scores, please install the following packages (for hg19):

``` r
BiocManager::install(c("BSgenome.Hsapiens.UCSC.hg19",
                       "fitCons.UCSC.hg19",
                       "phastCons100way.UCSC.hg19"))
```

Annotation of transcriptomic Features
=====================================

Introduction
------------

The ***WHISTLE*** R-package has been developed for the computational representation of the domain knowledge related to epitranscriptomic modifications, such as m6A, Pseudouridine, and m7G. The inputs of the main function ***transcriptomicFeatures*** are the `GRanges` object for range based RNA modification locations and the `TxDb` object for the transcript annotations with other optional annotation objects. The ***WHISTLE*** package has the following two key utilities:

-   Annotating the features that are powerful in the machine learning prediction of the mRNA modifications.

-   Elucidate the biological mechanisms of epitranscriptomic regulation by analyzing the statistical associations between the annotated features and other functional genomic features in field of epitranscriptomics.

Transcript annotation can be provided as a GTF/GFF file, a ***TxDb*** object, or automatically downloaded from UCSC through the internet using `makeTxDbFromUCSC`.

We will in the next see how the annotation of the transcriptomic features can be accomplished in a single command.

Generation of features on transcript topology
---------------------------------------------

Let us firstly load the WHISTLE package and the `TxDb` for transcript annotations.

``` r
library(WHISTLE)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
```

The main functioon of ***WHISTLE*** R-package is to annotate transcriptomic features based on a range based RNA modification sites input (Sites in single based resolution are highly recommended). The necessary inputs of this function are the RNA modification site in `GRanges` object and the transcript annotation files in `TxDb` object. In oder to extract features other than transcriptomic topologies, additional annotation objects such as `BSgenome` and `GSscores` are required. Please see `?transcriptomicFeatures` for detailed explainations.

``` r
feature_hg19 <- transcriptomicFeatures(
   mod_gr = m6A_ex_hg19,
   txdb = txdb_hg19
)
```

The output is a `data.frame` object with the columns indicating the annoated features in dummy variables or real valued number data types. The `data.frame` table can be used in prediction models directly as the feature/predictor matrix. The prediction models could be linear regression such as `lm` and `glm` in R base functions, or other advanced machine learning models such as svm and random forest implemented in R package `caret`.

``` r
str(feature_hg19) # Feature list as data.frame table
```

    ## 'data.frame':    100 obs. of  22 variables:
    ##  $ UTR5        : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##  $ UTR3        : logi  FALSE FALSE TRUE TRUE FALSE FALSE ...
    ##  $ stop_codon  : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##  $ start_codon : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##  $ TSS         : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##  $ exon_stop   : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ alt_exon    : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##  $ const_exon  : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
    ##  $ inter_exon  : logi  TRUE TRUE FALSE FALSE TRUE FALSE ...
    ##  $ long_exon   : logi  FALSE TRUE TRUE TRUE FALSE TRUE ...
    ##  $ last_ex400sc: logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##  $ pos_UTR5    : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ pos_UTR3    : num  0 0 0.685 0.319 0 ...
    ##  $ pos_exon    : num  0.201 0.445 0.69 0.799 0.217 ...
    ##  $ dist_5sj_2k : num  126 1325 2000 2000 93 ...
    ##  $ dist_3sj_2k : num  30 1060 2000 2000 24 277 8 98 19 192 ...
    ##  $ len_UTR3    : num  643 514 3882 504 0 ...
    ##  $ len_UTR5    : num  144 225 472 104 0 169 126 393 42 286 ...
    ##  $ len_gene_ex : num  1549 16254 7625 1706 9073 ...
    ##  $ ndist_m2h   : num  200 200 200 200 200 200 200 200 200 200 ...
    ##  $ sncRNA      : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...
    ##  $ lncRNA      : logi  FALSE FALSE FALSE FALSE FALSE FALSE ...

Generation of the complete feature set
--------------------------------------

In general, the features annotated by `transcriptomicFeatures` belong to the following categories:

-   ***Topological Features*** Features constructed from the transcriptomic landmarks. The toplogical features were previously reported by the published studies to associate with the locations of RNA modifiations.

-   ***Sequence Features*** The sequence derived features that are related to epitranscriptomics, such as the sequence motifs and the predicted RNA secondary structures.

-   ***Functional Genomic Features*** Features of other functional genomic information that have proofed correlation with the mRNA modifications, such as the functional types of the transcript, evolutionary conservation of the sequence, and the RBP binding sites of RNA modification writers, readers, and erasers.

Please notice that, when only the TxDb object is provided, only the topological features will be completely annotated in the returned `data.frame` table.

To generate the full feature set, please refer to the following code example for hg19.

``` r
library(BSgenome.Hsapiens.UCSC.hg19)
library(fitCons.UCSC.hg19)
library(phastCons100way.UCSC.hg19)

#Generate complete explanatory feature set for epitranscriptomic modifications

feature_hg19_complete <- transcriptomicFeatures(
   mod_gr = m6A_ex_hg19,
   
# TxDb object for the transcript annotation.
   txdb = txdb_hg19,

# BSgenome object for the genome sequence.
   bsgnm = Hsapiens,

# GScores object for the Fitness Consequence scores.
   fc = fitCons.UCSC.hg19,

# GScores object for the PhastCons scores.
   pc = phastCons100way.UCSC.hg19,

# RNA fold predicted RNA secondary structures for hg19 exons.
   struct_hybridize = Struc_hg19,

# Binding sites of m6A regulators.
   genomic_annot = RBP_m6A_hg19,

# Entrez gene ids of house keeping genes.
   hk_genes_id = hkGene_eid_hg19,

# Entrez gene ids of microRNA targeted genes.
   miRtar_genes_id = miRtarGene_eid_hg19,

# Instances of DRACH motif for m6A modification.
   motif = c("AAACA","AGACA","AAACT","AGACT","AAACC","AGACC",
             "GAACA","GGACA","GAACT","GGACT","GAACC","GGACC",
             "TAACA","TGACA","TAACT","TGACT","TAACC","TGACC"),

# The clustering effect on DRACH motifs.
   motif_clustering = "DRACH",

# The clustering on the input modification Granges.
   annot_clustering = m6A_ex_hg19
)
```

    ##    __  __  __ _______________    ______ 
    ##   / / / / / // / ___/_  __/ /   /  ___/ 
    ##  / / / / / // /\__ \ / / / /   /  ___/
    ## / /_/ /_/ // /___/ // / / /___/ /___/ 
    ## \___/\___//_/_____//_/ /_____/_____/     version 1.0.1 
    ## Type 'citation("WHISTLE")' for citing this R package in publications.
    ## Extracting exons, transcripts, and genes ... OK
    ## Handeling transcript isoform ambiguities ... OK
    ## Extracing the ambiguity removed exonic landmarkers ... OK
    ## Constructing transcriptomic features ... 
    ## Feature standardization: 'No Sandardization'.
    ## Isoform ambiguity method: 'Longest Mature RNA Transcript'.
    ## Gene ambiguity method: 'Average over Overlapping Genes'.
    ## The feature documentations are listed below:
    ## ============================================================================================================
    ##  ID  | Feature name |  Variable type   | Feature description 
    ## ============================================================================================================
    ##   1  |     UTR5     |  Dummy variable  | Site overlapping with 5'UTR.
    ## ------------------------------------------------------------------------------------------------------------
    ##   2  |     UTR3     |  Dummy variable  | Site overlapping with 3'UTR.
    ## ------------------------------------------------------------------------------------------------------------
    ##   3  |  stop_codon  |  Dummy variable  | Site within 201 bp centered by the stop codon.
    ## ------------------------------------------------------------------------------------------------------------
    ##   4  | start_codon  |  Dummy variable  | Site within 201 bp centered by the start codon.
    ## ------------------------------------------------------------------------------------------------------------
    ##   5  |     TSS      |  Dummy variable  | Site within 100bp downstream of TSS.
    ## ------------------------------------------------------------------------------------------------------------
    ##   6  |    TSS_A     |  Dummy variable  | Site within 100bp downstream of TSS that start with A.
    ## ------------------------------------------------------------------------------------------------------------
    ##   7  |  exon_stop   |  Dummy variable  | Site overlapping with exons containing stop codon.
    ## ------------------------------------------------------------------------------------------------------------
    ##   8  |   alt_exon   |  Dummy variable  | Site overlapping with alternatively spliced exons.
    ## ------------------------------------------------------------------------------------------------------------
    ##   9  |  const_exon  |  Dummy variable  | Site overlapping with constitutively spliced exons.
    ## ------------------------------------------------------------------------------------------------------------
    ##  10  |  inter_exon  |  Dummy variable  | Site overlapping with internal exons.
    ## ------------------------------------------------------------------------------------------------------------
    ##  11  |  long_exon   |  Dummy variable  | Site overlapping with long exons (exon length >= 400bp).
    ## ------------------------------------------------------------------------------------------------------------
    ##  12  | last_ex400sc |  Dummy variable  | Site within the 5' 400bp of the last exons containing stop codon.
    ## ------------------------------------------------------------------------------------------------------------
    ##  13  |   pos_UTR5   |  Real on [0,1]   | Relative positioning of the site on 5'UTR.
    ## ------------------------------------------------------------------------------------------------------------
    ##  14  |   pos_UTR3   |  Real on [0,1]   | Relative positioning of the site on 3'UTR.
    ## ------------------------------------------------------------------------------------------------------------
    ##  15  |   pos_exon   |  Real on [0,1]   | Relative positioning of the site on exon.
    ## ------------------------------------------------------------------------------------------------------------
    ##  16  | dist_5sj_2k  | Positive integer | Distance from site to the upstream (5' end) splicing junction.
    ## ------------------------------------------------------------------------------------------------------------
    ##  17  | dist_3sj_2k  | Positive integer | Distance from site to the downstream (3' end) splicing junction.
    ## ------------------------------------------------------------------------------------------------------------
    ##  18  |   len_UTR3   | Positive integer | Length of the 3'UTR overlapped by the site.
    ## ------------------------------------------------------------------------------------------------------------
    ##  19  |   len_UTR5   | Positive integer | Length of the 5'UTR overlapped by the site.
    ## ------------------------------------------------------------------------------------------------------------
    ##  20  | len_gene_ex  | Positive integer | length of the gene exonic regions overlapped by the site.
    ## ------------------------------------------------------------------------------------------------------------
    ##  21  |    AAACA     |  Dummy variable  | Site on sequence motif: AAACA.
    ## ------------------------------------------------------------------------------------------------------------
    ##  22  |    AGACA     |  Dummy variable  | Site on sequence motif: AGACA.
    ## ------------------------------------------------------------------------------------------------------------
    ##  23  |    AAACT     |  Dummy variable  | Site on sequence motif: AAACT.
    ## ------------------------------------------------------------------------------------------------------------
    ##  24  |    AGACT     |  Dummy variable  | Site on sequence motif: AGACT.
    ## ------------------------------------------------------------------------------------------------------------
    ##  25  |    AAACC     |  Dummy variable  | Site on sequence motif: AAACC.
    ## ------------------------------------------------------------------------------------------------------------
    ##  26  |    AGACC     |  Dummy variable  | Site on sequence motif: AGACC.
    ## ------------------------------------------------------------------------------------------------------------
    ##  27  |    GAACA     |  Dummy variable  | Site on sequence motif: GAACA.
    ## ------------------------------------------------------------------------------------------------------------
    ##  28  |    GGACA     |  Dummy variable  | Site on sequence motif: GGACA.
    ## ------------------------------------------------------------------------------------------------------------
    ##  29  |    GAACT     |  Dummy variable  | Site on sequence motif: GAACT.
    ## ------------------------------------------------------------------------------------------------------------
    ##  30  |    GGACT     |  Dummy variable  | Site on sequence motif: GGACT.
    ## ------------------------------------------------------------------------------------------------------------
    ##  31  |    GAACC     |  Dummy variable  | Site on sequence motif: GAACC.
    ## ------------------------------------------------------------------------------------------------------------
    ##  32  |    GGACC     |  Dummy variable  | Site on sequence motif: GGACC.
    ## ------------------------------------------------------------------------------------------------------------
    ##  33  |    TAACA     |  Dummy variable  | Site on sequence motif: TAACA.
    ## ------------------------------------------------------------------------------------------------------------
    ##  34  |    TGACA     |  Dummy variable  | Site on sequence motif: TGACA.
    ## ------------------------------------------------------------------------------------------------------------
    ##  35  |    TAACT     |  Dummy variable  | Site on sequence motif: TAACT.
    ## ------------------------------------------------------------------------------------------------------------
    ##  36  |    TGACT     |  Dummy variable  | Site on sequence motif: TGACT.
    ## ------------------------------------------------------------------------------------------------------------
    ##  37  |    TAACC     |  Dummy variable  | Site on sequence motif: TAACC.
    ## ------------------------------------------------------------------------------------------------------------
    ##  38  |    TGACC     |  Dummy variable  | Site on sequence motif: TGACC.
    ## ------------------------------------------------------------------------------------------------------------
    ##  39  |  ndist_m2h   | Positive integer | Distance to the nearest neighbooring sites (maximum to 200bp).
    ## ------------------------------------------------------------------------------------------------------------
    ##  40  |dist_DRACH_m2h| Positive integer | Distance to the nearest DRACH motif (maximum at 200bp).
    ## ------------------------------------------------------------------------------------------------------------
    ##  41  |    PC_1bp    |  Real on [0,1]   | Phast cons (PC) score of the site.
    ## ------------------------------------------------------------------------------------------------------------
    ##  42  |   PC_101bp   |  Real on [0,1]   | Average PC scores within 101bp window centered by the site.
    ## ------------------------------------------------------------------------------------------------------------
    ##  43  |    FC_1bp    |  Real on [0,1]   | Fitness consequence (FC) score of the site.
    ## ------------------------------------------------------------------------------------------------------------
    ##  44  |   FC_101bp   |  Real on [0,1]   | Average FC scores within 101bp window centered by the site.
    ## ------------------------------------------------------------------------------------------------------------
    ##  45  | struct_sterm |  Dummy variable  | Site within predicted sterm regions of RNA secondary structure.
    ## ------------------------------------------------------------------------------------------------------------
    ##  46  | struct_loop  |  Dummy variable  | Site within predicted loop regions of RNA secondary structure.
    ## ------------------------------------------------------------------------------------------------------------
    ##  47  |    METTL3    |  Dummy variable  | Site overlapping with range based annotation: METTL3.
    ## ------------------------------------------------------------------------------------------------------------
    ##  48  |   METTL14    |  Dummy variable  | Site overlapping with range based annotation: METTL14.
    ## ------------------------------------------------------------------------------------------------------------
    ##  49  |     WTAP     |  Dummy variable  | Site overlapping with range based annotation: WTAP.
    ## ------------------------------------------------------------------------------------------------------------
    ##  50  |    RBM15     |  Dummy variable  | Site overlapping with range based annotation: RBM15.
    ## ------------------------------------------------------------------------------------------------------------
    ##  51  |    RBM15B    |  Dummy variable  | Site overlapping with range based annotation: RBM15B.
    ## ------------------------------------------------------------------------------------------------------------
    ##  52  |    YTHDF1    |  Dummy variable  | Site overlapping with range based annotation: YTHDF1.
    ## ------------------------------------------------------------------------------------------------------------
    ##  53  |    YTHDF2    |  Dummy variable  | Site overlapping with range based annotation: YTHDF2.
    ## ------------------------------------------------------------------------------------------------------------
    ##  54  |    YTHDF3    |  Dummy variable  | Site overlapping with range based annotation: YTHDF3.
    ## ------------------------------------------------------------------------------------------------------------
    ##  55  |    YTHDC1    |  Dummy variable  | Site overlapping with range based annotation: YTHDC1.
    ## ------------------------------------------------------------------------------------------------------------
    ##  56  |    YTHDC2    |  Dummy variable  | Site overlapping with range based annotation: YTHDC2.
    ## ------------------------------------------------------------------------------------------------------------
    ##  57  |  HNRNPA2B1   |  Dummy variable  | Site overlapping with range based annotation: HNRNPA2B1.
    ## ------------------------------------------------------------------------------------------------------------
    ##  58  |    HNRNPC    |  Dummy variable  | Site overlapping with range based annotation: HNRNPC.
    ## ------------------------------------------------------------------------------------------------------------
    ##  59  |    ALKBH5    |  Dummy variable  | Site overlapping with range based annotation: ALKBH5.
    ## ------------------------------------------------------------------------------------------------------------
    ##  60  |     FTO      |  Dummy variable  | Site overlapping with range based annotation: FTO.
    ## ------------------------------------------------------------------------------------------------------------
    ##  61  |    sncRNA    |  Dummy variable  | Site overlapping with sncRNA (ncRNA <= 200bp).
    ## ------------------------------------------------------------------------------------------------------------
    ##  62  |    lncRNA    |  Dummy variable  | Site overlapping with lncRNA (ncRNA > 200bp).
    ## ------------------------------------------------------------------------------------------------------------
    ##  63  |   HK_genes   |  Dummy variable  | Site overlapping with house keeping genes.
    ## ------------------------------------------------------------------------------------------------------------
    ##  64  | miRtar_genes |  Dummy variable  | Site overlapping with microRNA targeted genes.
    ## ------------------------------------------------------------------------------------------------------------
    ## Transcriptomic features generation complete.

Contact
-------

Please contact the maintainer of WHISTLE if you have encountered any problems:

**ZhenWei**: <zhen.wei@xjtlu.edu.cn>

Session Info
------------

``` r
sessionInfo()
```

    ## R version 3.6.2 (2019-12-12)
    ## Platform: x86_64-apple-darwin15.6.0 (64-bit)
    ## Running under: macOS High Sierra 10.13.6
    ## 
    ## Matrix products: default
    ## BLAS:   /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRblas.0.dylib
    ## LAPACK: /Library/Frameworks/R.framework/Versions/3.6/Resources/lib/libRlapack.dylib
    ## 
    ## locale:
    ## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
    ## 
    ## attached base packages:
    ## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    ## [8] methods   base     
    ## 
    ## other attached packages:
    ##  [1] phastCons100way.UCSC.hg19_3.7.2        
    ##  [2] fitCons.UCSC.hg19_3.7.1                
    ##  [3] GenomicScores_1.10.0                   
    ##  [4] BSgenome.Hsapiens.UCSC.hg19_1.4.0      
    ##  [5] BSgenome_1.54.0                        
    ##  [6] rtracklayer_1.46.0                     
    ##  [7] Biostrings_2.54.0                      
    ##  [8] XVector_0.26.0                         
    ##  [9] TxDb.Hsapiens.UCSC.hg19.knownGene_3.2.2
    ## [10] GenomicFeatures_1.38.1                 
    ## [11] AnnotationDbi_1.48.0                   
    ## [12] Biobase_2.46.0                         
    ## [13] GenomicRanges_1.38.0                   
    ## [14] GenomeInfoDb_1.22.0                    
    ## [15] IRanges_2.20.2                         
    ## [16] S4Vectors_0.24.4                       
    ## [17] BiocGenerics_0.32.0                    
    ## [18] WHISTLE_0.99.1                        
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] httr_1.4.1                    bit64_0.9-7                  
    ##  [3] AnnotationHub_2.18.0          shiny_1.4.0                  
    ##  [5] assertthat_0.2.1              interactiveDisplayBase_1.24.0
    ##  [7] askpass_1.1                   BiocManager_1.30.10          
    ##  [9] BiocFileCache_1.10.2          blob_1.2.1                   
    ## [11] GenomeInfoDbData_1.2.2        Rsamtools_2.2.2              
    ## [13] yaml_2.2.1                    progress_1.2.2               
    ## [15] BiocVersion_3.10.1            pillar_1.4.3                 
    ## [17] RSQLite_2.2.0                 lattice_0.20-38              
    ## [19] glue_1.3.1                    digest_0.6.24                
    ## [21] promises_1.1.0                htmltools_0.4.0              
    ## [23] httpuv_1.5.2                  Matrix_1.2-18                
    ## [25] XML_3.99-0.3                  pkgconfig_2.0.3              
    ## [27] biomaRt_2.42.0                zlibbioc_1.32.0              
    ## [29] xtable_1.8-4                  purrr_0.3.3                  
    ## [31] later_1.0.0                   BiocParallel_1.20.1          
    ## [33] tibble_3.0.1                  openssl_1.4.1                
    ## [35] ellipsis_0.3.0                SummarizedExperiment_1.16.1  
    ## [37] magrittr_1.5                  crayon_1.3.4                 
    ## [39] mime_0.9                      memoise_1.1.0                
    ## [41] evaluate_0.14                 tools_3.6.2                  
    ## [43] prettyunits_1.1.1             hms_0.5.3                    
    ## [45] lifecycle_0.2.0               matrixStats_0.55.0           
    ## [47] stringr_1.4.0                 DelayedArray_0.12.3          
    ## [49] compiler_3.6.2                rlang_0.4.5                  
    ## [51] grid_3.6.2                    RCurl_1.98-1.1               
    ## [53] rappdirs_0.3.1                bitops_1.0-6                 
    ## [55] rmarkdown_2.1                 DBI_1.1.0                    
    ## [57] curl_4.3                      R6_2.4.1                     
    ## [59] GenomicAlignments_1.22.1      knitr_1.28                   
    ## [61] dplyr_0.8.4                   fastmap_1.0.1                
    ## [63] bit_1.1-15.2                  stringi_1.4.5                
    ## [65] Rcpp_1.0.3                    vctrs_0.2.4                  
    ## [67] dbplyr_1.4.2                  tidyselect_1.0.0             
    ## [69] xfun_0.12
