#' @title Generation of transcriptomic features for range based epitranscriptomic dataset
#'
#' @description The function \bold{\code{transcriptomicFeatures}} could generate features based on the \code{\link{GRanges}} input object of the sites for RNA epitranscriptomic modification.
#' Users could provide additional genomic and transcripomic annotations to expand the feature sets annotated by this function.
#'
#' Under the default setting, \code{transcriptomicFeatures} could provide 40 topological features based on the \code{\link{TxDb}} object of the transcript annotation.
#' More than 80 features related to the domain knowledge of epitranscriptomic modifications could be annotated given a complete support of the annotation objects. The suggested additonal annotation data types include:
#'
#' \enumerate{
#' \item \code{\link{BSgenome}} object for the genome sequence information, such as the \code{BSgenome.Hsapiens.UCSC.hg19} package that provides the genomic sequences for the hg19 genome assembly.
#'
#' \item \code{\link{GScores}} objects for the annotation of evolutionary conservation scores, such as the \code{phastCons100way.UCSC.hg19} package that provides the annotation of Phastcons scores for hg19.
#'
#' \item \code{GRanges} or \code{\link{GRangesList}} objects that providing the annotations of other range based transcriptomic features gathered from diffrent databases.
#' }
#'
#' Please download the suggested annotation packages from \href{https://www.bioconductor.org/}{Bioconductor}.
#'
#' @param mod_gr A \code{GRanges} object for the genomic location of RNA modification sites.
#'
#' @param txdb A \code{TxDb} object for the transcript annotation.
#'
#' The \code{TxDb} can be obtained from either Bioconductor or from the GFF/GTF files using the function \code{\link{makeTxDbFromGFF}}.
#'
#' @param bsgnm Optional; a \code{BSgenome} object for the genomic sequence annotation.
#'
#' The \code{BSgenome} can be downloaded either from Bioconductor or through the function \code{\link{getBSgenome}}.
#'
#' @param fc,pc Optional; \code{GScores} objects for the annotation of phastCons[1] scores and the fitness consequences scores[2].
#'
#' The \code{GSscores} can be downloaded from Bioconductor.
#'
#' @param struct_hybridize A \code{GRanges} object recording the hybridized region of the RNA transcripts.
#'
#' WHISTLE includes the predicted exonic RNA secondary structures of hg19 and mm10 in the exported data.
#'Please check \code{\link{Struc_hg19}} and \code{\link{Struc_mm10}} for details.
#'
#' @param genomic_annot A \code{list} of \code{\link{GRanges}} for the user defined range based functional genomic annotation.
#'
#' New dummy variable features will be created for each list element by overlapping the \code{GRanges} with the corresponding range based annotations.
#' The names of the list will correspond to the names of features in the returned \code{data.frame} object.
#'
#' Some of the range based functional genomic annotations can be important in the site prediction of RNA modifications.
#' These annotations usually include the binding sites of the regulator proteins (e.x. writers, erasers, readers) related to the RNA modification.
#' Additionally, the microRNA binding targets can associate with the locations of certain RNA modification.
#'
#' WHISTLE package includes the following build-in \code{list} objects to represent the RBP binding sites of the regulator proteins for different RNA modifications:
#' \code{\link{RBP_m6A_hg19}}, \code{\link{RBP_m6A_mm10}}, \code{\link{RBP_m1A_hg19}}, \code{\link{RBP_m5C_hg19}}, \code{\link{RBP_Psi_hg19}}; all of the RBP data in the list objects are gathered from Starbase[5].
#'
#' Furthermore, the information of microRNA binding sites could be found in :\code{\link{TargetScan_hg19}}.
#'
#' @param motif A character vector for the sequence motifs centered by the RNA modification site.
#' Each element of the character vector will become an annotated feature in the returned table.
#'
#' For example, to predict m6A RNA modification, the motif vector could be the combination of the instances in RRACH.
#'
#' P.S. The motif will not be attached if the \code{mod_gr} is not in single nucleotide resolution (with all width = 1).
#'
#' @param annot_clustering A \code{GRanges} object to generate clustering indexes of the modification sites related to a particular range based annotation;
#' under the default setting, the clustering will be conducted on the input GRanges.
#'
#' 4 features --- \code{clust_f100}, \code{clust_f1000}, \code{dist_nearest_p200}, and \code{dist_nearest_p2000} will be generated to represent the clustering effect of the provided annotation related to the RNA modification site.
#'
#' @param motif_clustering A \code{character} indicating the motif used when generating the features for the motif clustering; default: "DRACH".
#'
#' The motif clustering will represent the clustering of the sequence motif using 4 features similar to the annotation clusering.
#'
#' @param hk_genes_id Optional; A character string for the gene id of the house keeping genes. The ids should correspond to the gene ids in the provided \code{TxDb} object.
#'
#' The entrez gene IDs of the house keeping genes for hg19 are included in WHISTLE: \code{\link{HK_hg19_eids}} [6].
#'
#' @param miRtar_genes_id  Optional; A character string for the Gene IDs of the microRNA targeted genes. The id should also correspond to the ids in the provided \code{TxDb} object.
#'
#' The entrez gene IDs of the microRNA targeted genes for hg19 are included in WHISTLE: \code{\link{miRtarGene_eid_hg19}}.
#'
#' @param isoform_ambiguity_method Can be \code{"longest_tx"} or \code{"average"}. The former keeps only the longest transcript as the transcript annotation.
#' The later will use the average feature entries for multiple mappings of the transcript isoform; default: \code{"longest_tx"}.
#'
#' @param gene_ambiguity_method Can be \code{"drop_overlap"} or \code{"average"}. The former will not annotate the modification Site overlapped with > 1 genes (By returning NA).
#' The later will use the average feature entries for mapping of multiple genes; default: \code{"average"}.
#'
#' @param standardization A logical indicating whether to standardize the continous features. For a given scalar feature X, the standardization is conducted by \eqn{ \frac{X - mean(X)}{sd(X)} }{(X - mean(X))/sd(X)}; default: \code{FALSE}.
#'
#' @details \code{transcriptomicFeatures} in \bold{WHISTLE} package can retrieve a set of domain knowledge based transcriptomic features that contribute toward the topology of epitranscriptomic modifications.
#'
#' The utility of the annotated transcriptomic features is 2 folds:
#'
#' \enumerate{
#' \item The transcriptomic features are useful in the machine learning prediction of the epitranscriptomic markers.
#' WHISTLE predictor utilizes part of the transcriptomic features to achieve superior performances on the m6A site preidiction than any other predictors based only on the sequence features[3].
#'
#' \item The features are also important for the topological characterization of the epitranscriptomic regulations under a new biological condition.
#' The statistical association (e.x. with linear or logistic regression) between the transcriptomic feature and the RNA modification measurements can elucidate the scientific significance of a particular feature in epitranscriptomic regulations[4].
#'
#' }
#'
#' In general, \code{transcriptomicFeatures} could annotate the transcriptomic features that belong to the following categories:
#'
#' \describe{
#' \item{\strong{Topological features}}{Features constructed from the transcriptomic landmarks. The toplogical features were previously reported by the published studies to associate with the locations of RNA modifiations.}
#'
#' \item{\strong{Sequence features}}{The sequence derived features that are related to epitranscriptomics, such as the sequence motifs, predicted RNA secondary structures, and GC contents.}
#'
#' \item{\strong{Functional genomic features}}{Features of other functional genomic information that have proofed correlation with the mRNA modifications,
#' such as the functional aspects of the transcript, evolutionary conservation of the sequence, and the RBP binding sites of RNA modification writers, readers, and erasers. }
#' }
#'
#' By default, only the topological features will be annotated if users provide the \code{TxDb} object along without additional annotation objects.
#'
#' @return  The function accepts the input of a \code{GRanges} object for RNA modification sites, and it will return a \code{data.frame} with columns of the annotated transcriptomic features.
#' The rows of the returned \code{data.frame} object will correspond to the order of sites in the corresponding GRanges.
#'
#' @examples
#'
#' ## --------------------------------------------------------------------
#' ## Topological features for hg19
#' ## --------------------------------------------------------------------
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#'
#' txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' feature_hg19 <- transcriptomicFeatures(
#' mod_gr = m6A_ex_hg19,
#' txdb = txdb_hg19
#' )
#'
#' str(feature_hg19)
#'
#' ## --------------------------------------------------------------------
#' ## Topological features for mm10
#' ## --------------------------------------------------------------------
#' \dontrun{
#' library(TxDb.Mmusculus.UCSC.mm10.knownGene)
#'
#' txdb_mm10 <- TxDb.Mmusculus.UCSC.mm10.knownGene
#'
#' feature_mm10 <- transcriptomicFeatures(
#' mod_gr = m6A_ex_mm10,
#' txdb = txdb_mm10
#' )
#'
#' str(feature_mm10)
#'
#' ## --------------------------------------------------------------------
#' ## Complete features for hg19
#' ## --------------------------------------------------------------------
#'
#' library(TxDb.Hsapiens.UCSC.hg19.knownGene)
#' library(BSgenome.Hsapiens.UCSC.hg19)
#' library(fitCons.UCSC.hg19)
#' library(phastCons100way.UCSC.hg19)
#'
#' txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
#'
#' feature_complete <- transcriptomicFeatures(
#'
#'    mod_gr = m6A_ex_hg19,
#
#' # TxDb object for the transcript annotation.
#'    txdb = txdb_hg19,
#'
#' # BSgenome object for the genome sequence.
#'    bsgnm = Hsapiens,
#'
#' # GScores object for the Fitness Consequence scores.
#'    fc = fitCons.UCSC.hg19,
#'
#' # GScores object for the PhastCons scores.
#'    pc = phastCons100way.UCSC.hg19,
#'
#' # RNA fold predicted RNA secondary structures for hg19 exons.
#'    struct_hybridize = Struc_hg19,
#'
#' # Binding sites of m6A regulators.
#'    genomic_annot = RBP_m6A_hg19,
#'
#' # Entrez gene ids of house keeping genes.
#'    hk_genes_id = hkGene_eid_hg19,
#'
#' # Entrez gene ids of microRNA targeted genes.
#'    miRtar_genes_id = miRtarGene_eid_hg19,
#'
#' # Instances of DRACH motif for m6A modification.
#'    motif = c("AAACA","AGACA","AAACT","AGACT","AAACC","AGACC",
#'              "GAACA","GGACA","GAACT","GGACT","GAACC","GGACC",
#'              "TAACA","TGACA","TAACT","TGACT","TAACC","TGACC"),
#'
#' # The clustering effect on DRACH motifs.
#'    motif_clustering = "DRACH",
#'
#' # The clustering on the input modification Granges.
#'    annot_clustering = m6A_ex_hg19
#' )
#'
#' str(feature_complete)
#'
#' }
#' @references
#'
#' \enumerate{
#'
#' \item Siepel A, Bejerano G, Pedersen JS, Hinrichs AS, Hou M, Rosenbloom K, Clawson H, Spieth J, Hillier LW, Richards S: Evolutionarily conserved elements in vertebrate, insect, worm, and yeast genomes. Genome research 2005, 15(8):1034-1050.
#'
#' \item Gulko B, Hubisz MJ, Gronau I, Siepel A: A method for calculating probabilities of fitness consequences for point mutations across the human genome. Nature genetics 2015, 47(3):276.
#'
#' \item Chen K, Wei Z, Zhang Q, Wu X, Rong R, Lu Z, Su J, de Magalhaes JP, Rigden DJ, Meng J: WHISTLE: a high-accuracy map of the human N6-methyladenosine (m6A) epitranscriptome predicted using a machine learning approach. Nucleic Acids Research 2019.
#'
#' \item Schwartz S, Mumbach M, Jovanovic M, Wang T, Maciag K, Bushkin GG, Mertins P, Ter-Ovanesyan D, Habib N, Cacchiarelli D: Perturbation of m6A Writers Reveals Two Distinct Classes of mRNA Methylation at Internal and 5' Sites. Cell Reports 2014, 8(1):284-296.
#'
#' \item Li J-H, Liu S, Zhou H, Qu L-H, Yang J-H: starBase v2.0: decoding miRNA-ceRNA, miRNA-ncRNA and protein–RNA interaction networks from large-scale CLIP-Seq data. Nucleic acids research 2013, 42(D1):D92-D97.
#'
#' \item Eisenberg E, Levanon EY: Human housekeeping genes, revisited. Trends in Genetics 2013, 29(10):569-574.
#'
#' }
#'
#'
#' @import BSgenome
#' @import Biostrings
#' @import GenomicFeatures
#' @import GenomicRanges
#' @importFrom GenomicScores score
#' @export
#' 
transcriptomicFeatures <- function(mod_gr,
                            txdb,
                            bsgnm = NULL,
                            fc = NULL,
                            pc = NULL,
                            struct_hybridize = NULL,
                            genomic_annot = NULL,
                            motif = NULL,
                            annot_clustering = NULL,
                            motif_clustering = NULL,
                            hk_genes_id = NULL,
                            miRtar_genes_id = NULL,
                            isoform_ambiguity_method = c("longest_tx",
                                                          "average"),
                            gene_ambiguity_method = c("average",
                                                       "drop_overlap"),
                            standardization = FALSE
) {

  #Pre_check
  stopifnot(length(mod_gr) != 0)
  isoform_ambiguity_method <- match.arg(isoform_ambiguity_method)
  gene_ambiguity_method <- match.arg(gene_ambiguity_method)

  cat("   __  __  __ _______________    ______ \n")
  cat("  / / / / / // / ___/_  __/ /   /  ___/ \n" )
  cat(" / / / / / // /\\__ \\ / / / /   /  ___/ \n")
  cat("/ /_/ /_/ // /___/ // / / /___/ /___ \n")
  cat("\\___/\\___//_/_____//_/ /_____/_____/   version 1.0.1 \n")
  cat("Type 'citation(\"WHISTLE\")' for citing this R package in publications.\n")


  #Extract genomic features.
  cat("Extracting exons, transcripts, and genes ... ",append = T)
  exbytx_txdb <- exonsBy(txdb, by = "tx")
  tx <- transcripts(txdb)
  names(tx) <- tx$tx_id
  cat("OK\n")

  #Handel the isoform ambiguity.

  if (isoform_ambiguity_method == "longest_tx") {

    cat("Handeling transcript isoform ambiguities ... ",append = T)
    Longest_tx_table <- find_longest_transcript(exbytx_txdb, txdb)
    Kept_tx_indx <- Longest_tx_table$TXID[Longest_tx_table$longest]
    rm(Longest_tx_table)
    cat("OK\n")

  } else {
    Kept_tx_indx <- T
  }

  #Filter the transcripts into the longest ones by their genes
  exbytx_txdb <- exbytx_txdb[Kept_tx_indx]
  tx <- tx[Kept_tx_indx]

  #Remove the overlapped transcripts that belong to multiple genes.
  if (gene_ambiguity_method == "drop_overlap") {
    cat("Removing overlapped transcripts that belong to multiple genes ... ",append = T)
    exbytx_txdb <-
      exbytx_txdb[countOverlaps(exbytx_txdb, exbytx_txdb) == 1]
    cat("OK\n")
  }

  #PS. exbg_txdb is used independently of the isoform ambiguity (but the gene dis-ambiguity is still used).

  cat("Extracing the ambiguity removed exonic landmarkers ... ",append = T)

  exbg_txdb <- exonsBy(txdb, by = "gene")

  if (gene_ambiguity_method == "drop_overlap") {
    keep_indx <- countOverlaps(exbg_txdb, exbg_txdb) == 1
    removed_exbg_txdb <- exbg_txdb[!keep_indx]
    exbg_txdb <- exbg_txdb[keep_indx]
    rm(keep_indx)
  }

  #Extract all the ambiguity removed exonic regions.

  exs_txdb <- unlist(exbytx_txdb)

  txid <- names(exbytx_txdb)

  cds <- cdsBy(txdb, by = "tx")

  cds <- cds[names(cds) %in% txid]

  Stop_codons <- resize(unlist(range(cds)), 3, fix = "end")

  Start_codons <- resize(unlist(range(cds)), 3, fix = "start")

  TSS <- resize(unlist(range(exbytx_txdb)), 1, fix = "start")

  TSS <- resize(TSS, 101, fix = "start")

  UTR3 <- threeUTRsByTranscript(txdb)

  UTR3 <- UTR3[names(UTR3) %in% txid]

  UTR5 <- fiveUTRsByTranscript(txdb)

  UTR5 <- UTR5[names(UTR5) %in% txid]

  Feature_matrix = data.frame(UTR5 = mod_gr %over% UTR5)

  Intron <- intronsByTranscript(txdb)

  #Define a function for standardizing features
  scale_i <- function(x, stand = T) {
    if (stand)
      return((x - mean(x)) / sd(x))
    else
      return(x)
  }
  cat("OK\n")

  #Report Features

  cat("Constructing transcriptomic features ... \n")


  cat(paste0("Feature standardization: '", ifelse(standardization,"Z Scores for Numeric Features","No Sandardization"), "'.\n"))

  cat(paste0("Isoform ambiguity method: ", ifelse(isoform_ambiguity_method == "longest_tx",
                                                  "'Longest Mature RNA Transcript'",
                                                  "'Average over Isoforms'"), ".\n"))

  cat(paste0("Gene ambiguity method: ", ifelse(gene_ambiguity_method == "drop_overlap",
                                                "'Drop Overlapping Genes'",
                                                "'Average over Overlapping Genes'"),".\n"))


  cat("The feature documentations are listed below:\n")

  cat("============================================================================================================\n")
  cat(" ID  | Feature name |  Variable type   | Feature description \n")
  cat("============================================================================================================\n")

  i = 1

  ## Define the helper functions reporting the progress

  render_table <- function(char, expected_length){
   round_length <- (expected_length - nchar(char))/2
   paste0( paste( rep(" ",floor(round_length)+1), collapse = ""),char,
               paste( rep(" ",ceiling(round_length)+1), collapse = ""),"|")
  }

  Speak <-
    function( I, variable_name, variable_type, feature_detail, last = F) {
      cat(paste0(render_table(I, 3),
                 render_table(variable_name, 12),
                 render_table(variable_type, 16),
                 " ", feature_detail, "\n"))
      if(last){
        cat("============================================================================================================\n")
      }else{
        cat("------------------------------------------------------------------------------------------------------------\n")
      }
      return(I + 1)
    }

  i  = Speak(i,"UTR5","Dummy variable","Site overlapping with 5'UTR.")

  Feature_matrix$UTR3 <- mod_gr %over% UTR3

  i  = Speak(i, "UTR3", "Dummy variable", "Site overlapping with 3'UTR.")

  # - Stop_codons: Stop codon (201 bp center).

  Feature_matrix$stop_codon <- mod_gr %over% (Stop_codons + 100)

  i  = Speak(i, "stop_codon", "Dummy variable", "Site within 201 bp centered by the stop codon.")

  # - Start_codons: Start codon (201 bp center).
  #

  Feature_matrix$start_codon <- mod_gr %over% (Start_codons + 100)
  i  = Speak(i, "start_codon", "Dummy variable", "Site within 201 bp centered by the start codon.")

  # - m6Am: 5'Cap m6Am (TSS that has underlying sequence of A).
  #
  Feature_matrix$TSS <- mod_gr %over% TSS

  i  = Speak(i, "TSS", "Dummy variable", "Site within 100bp downstream of TSS.")

  # - m6Am: 5'Cap m6Am (TSS that has underlying sequence of A).
  #
  if(!is.null(bsgnm)){
    A_idx <- vcountPattern("A", DNAStringSet(Views(bsgnm, TSS))) > 0
    Feature_matrix$TSS_A <- mod_gr %over% TSS[A_idx]
    rm(A_idx)
    i = Speak(i, "TSS_A", "Dummy variable", "Site within 100bp downstream of TSS that start with A.")
  }

  # - Annotate various exonic regions.
  Feature_matrix$exon_stop <-
    mod_gr %over% subsetByOverlaps(exbg_txdb,
                                   Stop_codons)

  i  = Speak(i, "exon_stop", "Dummy variable", "Site overlapping with exons containing stop codon.")

  #Alternative exons: exons that could be introns in some transcripts
  exbytx_unlist <- unlist(exbytx_txdb)

  Feature_matrix$alt_exon <-
    mod_gr %over% subsetByOverlaps(exs_txdb,
                                   unlist(Intron) + 1,
                                   type = "within",
                                   maxgap =
                                     0L)

  i  = Speak(i, "alt_exon", "Dummy variable", "Site overlapping with alternatively spliced exons.")

  #Constitutive exons: exons that always present in all transcripts

  Feature_matrix$const_exon <-
    mod_gr %over% subsetByOverlaps(
      exbytx_unlist,
      unlist(Intron) + 1,
      type = "within",
      maxgap = 0L,
      invert = T
    )

  i  = Speak(i, "const_exon", "Dummy variable", "Site overlapping with constitutively spliced exons.")

  #Internal exons: exons that are not the first or the last exons
  ex_ranks <- exbytx_unlist$exon_rank
  names(ex_ranks) = 1:length(ex_ranks)

  Idx_last_exon <-
    tapply(ex_ranks , names(exbytx_unlist) , function(x)
      names(x)[which.max(x)])
  Idx_last_exon <-
    as.numeric(Idx_last_exon[unique(names(exbytx_unlist))])
  Indx_last_exon <- vector("logical", length(ex_ranks))
  Indx_last_exon[Idx_last_exon] <- T

  Indx_first_exon <- ex_ranks == 1
  rm(ex_ranks)

  Feature_matrix$inter_exon <-
    mod_gr %over% exbytx_unlist[!(Indx_last_exon | Indx_first_exon)]

  rm(Indx_first_exon)

  i  = Speak(i, "inter_exon", "Dummy variable", "Site overlapping with internal exons.")

  #
  # - long_exons: Long exon (length > 400bp).
  #

  long_exs_txdb <- exs_txdb[width(exs_txdb) >= 400]
  Feature_matrix$long_exon <- mod_gr %over% long_exs_txdb
  rm(long_exs_txdb)

  i  = Speak(i, "long_exon", "Dummy variable", "Site overlapping with long exons (exon length >= 400bp).")

  #
  ### below tries to annotate the features proposed by Ke's paper
  #

  last_exons <- exbytx_unlist[Indx_last_exon]

  rm(Indx_last_exon)

  last_exons_400bp <-
    resize(last_exons[width(last_exons) >= 400], 400, fix = "start")

  Feature_matrix$last_ex400sc <-
    mod_gr %over% subsetByOverlaps(last_exons_400bp, Stop_codons)

  rm(Stop_codons, last_exons_400bp)

  i  = Speak(i, "last_ex400sc", "Dummy variable", "Site within the 5' 400bp of the last exons containing stop codon.")

  #
  #
  #  ###2. Relative positions### |- The entries fall into the scale of [0,1], if the entries not mapped to the range, the value is 0.
  #
  # - pos_UTR5: Relative positioning on 5'UTR.

  Feature_matrix$pos_UTR5 <-
    relative_pos_map(mod_gr, UTR5, 0, standardization)

  i  = Speak(i, "pos_UTR5", "Real on [0,1]", "Relative positioning of the site on 5'UTR.")

  #
  # - pos_UTR3: Relative positioning on 3'UTR.
  #

  Feature_matrix$pos_UTR3 <-
    relative_pos_map(mod_gr, UTR3, 0, standardization)

  i  = Speak(i, "pos_UTR3", "Real on [0,1]", "Relative positioning of the site on 3'UTR.")

  # - pos_exons: Relative positioning on exons.

  Feature_matrix$pos_exon <-
    relative_pos_map(mod_gr, exs_txdb, 0, standardization)

  i  = Speak(i, "pos_exon", "Real on [0,1]", "Relative positioning of the site on exon.")


  #
  # ### Splicing. ###
  #

  Splicing_junctions <-
    reduce(c(
      resize(exbytx_unlist, 1, fix = "start"),
      resize(exbytx_unlist, 1, fix = "end")
    ), min.gapwidth = 0L)

  Splicing_junctions <-
    subsetByOverlaps(Splicing_junctions , c(resize(TSS, 1, fix = "start"),
                                            resize(unlist(range(
                                              exbytx_txdb
                                            )), 1, fix = "end")), invert = T)

  rm(TSS)

  Feature_matrix$dist_5sj_2k <- distance_map(
    mod_gr,
    Splicing_junctions,
    end = "five",
    maximum_distance = 2000,
    standardize = standardization
  )

  i  = Speak(i, "dist_5sj_2k", "Positive integer", "Distance from site to the upstream (5' end) splicing junction.")

  Feature_matrix$dist_3sj_2k <- distance_map(
    mod_gr,
    Splicing_junctions,
    end = "three",
    maximum_distance = 2000,
    standardize = standardization
  )

  rm(Splicing_junctions)

  i  = Speak(i, "dist_3sj_2k", "Positive integer", "Distance from site to the downstream (3' end) splicing junction.")

  #
  #   ###3. Region length###
  #
  # - long_UTR3: Long 3'UTR (length > 400bp).


  Feature_matrix$len_UTR3 <- properties_map(
    query_gr = mod_gr,
    feature_gr = tx[names(UTR3)],
    feature_properties = sum(width(UTR3)),
    no_map_val = NA,
    normalize = standardization
  )

  rm(UTR3)

  Feature_matrix$len_UTR3[is.na(Feature_matrix$len_UTR3)] = 0

  i  = Speak(i, "len_UTR3", "Positive integer", "Length of the 3'UTR overlapped by the site.")

  Feature_matrix$len_UTR5 <- properties_map(
    query_gr = mod_gr,
    feature_gr = tx[names(UTR5)],
    feature_properties = sum(width(UTR5)),
    no_map_val = NA,
    normalize = standardization
  )

  rm(UTR5)

  Feature_matrix$len_UTR5[is.na(Feature_matrix$len_UTR5)] = 0

  i  = Speak(i, "len_UTR5", "Positive integer", "Length of the 5'UTR overlapped by the site.")

  # - Gene_length_ex: standardized gene length of exonic regions (z score).
  #

  txbg_txdb <- transcriptsBy(txdb, by = "gene")

  Feature_matrix$len_gene_ex <- properties_map(
    query_gr = mod_gr,
    feature_gr = txbg_txdb,
    feature_properties = sum(width(reduce(exbg_txdb))),
    no_map_val = NA,
    normalize = standardization
  )

  Feature_matrix$len_gene_ex[is.na(Feature_matrix$len_gene_ex)] = 0

  i  = Speak(i, "len_gene_ex", "Positive integer", "length of the gene exonic regions overlapped by the site.")
  rm(tx)

  #  #####=============== The following features that are optional ===============#####
  #
  #   ###4. Motif###
  #
  #
  if (is.null(motif)){ }else{
  if (is.null(bsgnm) & !is.null(motif)) {
    warning("BSgenome object not provided, the motif related features are not attached.")
  } else {
  if (any(width(mod_gr) != 1)) {
    warning("At least 1 range with width > 1, the motifs are not attached.")
  } else{
    is_motif <-
      function(motif, dss, exact = F)
        vcountPattern(DNAString(motif), dss, fixed = exact) > 0
    DSS <- DNAStringSet(Views(bsgnm, mod_gr + 2))

    for (m_i in motif)  {
      Feature_matrix[[m_i]] <- is_motif(m_i, DSS, F)
      i  = Speak(i, m_i, "Dummy variable", paste0("Site on sequence motif: ",m_i, "."))
    }
  }
  }
  }

  #
  # ## Clustering effect.
  #
  if (!is.null(annot_clustering)) {

    dist_self <- rep(2000, length(mod_gr))

    match_obj <- distanceToNearest(mod_gr, annot_clustering)

    dist_self[queryHits(match_obj)] <- mcols(match_obj)$distance

    rm(match_obj)

    #For mod_gr that is overlapping  with annot_clustering
    if (any(mod_gr %over% annot_clustering)) {
      match_obj <- distanceToNearest(annot_clustering)

      fol <- findOverlaps(mod_gr, annot_clustering)

      dist_tmp <- rep(2000, queryLength(fol))

      dist_tmp <-
        mcols(match_obj)$distance[match(subjectHits(fol), queryHits(match_obj))]

      dist_self[queryHits(fol)] <- dist_tmp

      dist_self[is.na(dist_self)] <- median(dist_self,na.rm = T)

      rm(fol, match_obj, dist_tmp)

    }

    Feature_matrix$ndist_m2h <-
      as.numeric(scale_i(pmin(dist_self, 200) , standardization))

    i  = Speak(i, "ndist_m2h", "Positive integer",  "Distance to the nearest neighbooring sites (maximum to 200bp).")

    rm(dist_self)

  }else{}

  if (!is.null(motif_clustering) & !is.null(bsgnm)) {
    tx_reduced <- reduce(transcripts(txdb))

    mod_gr_expanded <- reduce(mod_gr + 2000)

    motif_regions <- intersect(mod_gr_expanded, tx_reduced)

    rm(mod_gr_expanded, tx_reduced)

    motif_transcripts <- sample_sequence(motif_clustering,
                                         motif_regions,
                                         bsgnm)

    rm(motif_regions)


    motif_transcripts <-
      subsetByOverlaps(motif_transcripts, mod_gr, invert = T) #Remove all the self motifs

    match_obj <- distanceToNearest(mod_gr, motif_transcripts)

    dist_motif <- rep(2000, length(mod_gr))

    dist_motif[queryHits(match_obj)] <- mcols(match_obj)$distance

    dist_self <- rep(2000, length(mod_gr))

    match_obj <- distanceToNearest(mod_gr)

    dist_self[queryHits(match_obj)] <- mcols(match_obj)$distance

    less_index <- dist_self < dist_motif

    dist_motif[less_index] <- dist_self[less_index]


    Feature_matrix[[paste0("dist_", motif_clustering, "_m2h")]] <-
      as.numeric(scale_i(pmin(dist_motif, 200) , standardization))

    i  = Speak(i, paste0("dist_", motif_clustering, "_m2h"),
               "Positive integer",
               paste0(
                 "Distance to the nearest ",
                 motif_clustering,
                 " motif (maximum at 200bp)."
               ))

    rm(match_obj)
    rm(dist_self)
    rm(dist_motif)
    rm(less_index)
    rm(motif_transcripts)
  }

  #   ###5. Evolutionary fitness###
  #
  # - PC 1bp: standardized PC score 1 nt.

  if (!is.null(pc)) {
    suppressMessages(suppressWarnings(Feature_matrix$PC_1bp <-
                                        score(pc, mod_gr)))

    Feature_matrix$PC_1bp[is.na(Feature_matrix$PC_1bp)] = mean(na.omit(Feature_matrix$PC_1bp))

    Feature_matrix$PC_1bp = as.numeric(scale_i(Feature_matrix$PC_1bp , standardization))

    i  = Speak(i, paste0("PC_1bp"),
               "Real on [0,1]",
               "Phast cons (PC) score of the site.")

    suppressMessages(suppressWarnings(Feature_matrix$PC_101bp <-
                                        score(pc, mod_gr + 50)))

    Feature_matrix$PC_101bp[is.na(Feature_matrix$PC_101bp)] = mean(na.omit(Feature_matrix$PC_101bp))

    Feature_matrix$PC_101bp = as.numeric(scale_i(Feature_matrix$PC_101bp , standardization))

    i  = Speak(i, paste0("PC_101bp"),
               "Real on [0,1]",
               "Average PC scores within 101bp window centered by the site.")
  }

  #Feature 21. fitness consequences scores 1bp.

  if (!is.null(fc)) {
    suppressMessages(suppressWarnings(Feature_matrix$FC_1bp <-
                                        score(fc, mod_gr)))

    Feature_matrix$FC_1bp[is.na(Feature_matrix$FC_1bp)] = mean(na.omit(Feature_matrix$FC_1bp))

    Feature_matrix$FC_1bp = as.numeric(scale_i(Feature_matrix$FC_1bp , standardization))

    i  = Speak(i, paste0("FC_1bp"),
               "Real on [0,1]",
               "Fitness consequence (FC) score of the site.")

    suppressMessages(suppressWarnings(Feature_matrix$FC_101bp <-
                                        score(fc, mod_gr + 50)))

    Feature_matrix$FC_101bp[is.na(Feature_matrix$FC_101bp)] = mean(na.omit(Feature_matrix$FC_101bp))

    Feature_matrix$FC_101bp =  as.numeric(scale_i(Feature_matrix$FC_101bp , standardization))

    i  = Speak(i, paste0("FC_101bp"),
               "Real on [0,1]",
               "Average FC scores within 101bp window centered by the site.")
  }

  if (is.null(struct_hybridize)) {
  } else {
    Feature_matrix$struct_sterm <- mod_gr %over% struct_hybridize

    i  = Speak(i, "struct_sterm",
               "Dummy variable",
               "Site within predicted sterm regions of RNA secondary structure.")

    Feature_matrix$struct_loop <-
      mod_gr %over% infer_loop(struct_hybridize)

    i  = Speak(i, "struct_loop",
               "Dummy variable",
               "Site within predicted loop regions of RNA secondary structure.")

  }

  #
  #   ###6. User specified features by argument \code{genomic_annot}###
  #
  # The entries are logical / dummy variables, specifying whether overlapping with each GRanges or GRanges list.

  if (!is.null(genomic_annot)) {
    for (i_flst in names(genomic_annot)) {
      suppressWarnings(Feature_matrix[[i_flst]] <-
                         mod_gr %over% genomic_annot[[i_flst]])
      i  = Speak(i, i_flst,
                 "Dummy variable",
                 paste0("Site overlapping with range based annotation: ", i_flst, "."))

      }
  }

  #
  #   ###7.Gene attribute###


  #
  # - sncRNA: small noncoding RNA (<= 200bp)
  #

  coding_tx <- names(cds)
  exbytx_nc <- exbytx_txdb[!names(exbytx_txdb) %in% coding_tx]
  lnc_idx <- sum(width(exbytx_nc)) > 200
  Feature_matrix$sncRNA  <- mod_gr %over% exbytx_nc[!lnc_idx]
  rm(coding_tx)

  i  = Speak(i, "sncRNA",
             "Dummy variable",
             paste0("Site overlapping with sncRNA (ncRNA <= 200bp)."))

  # - lncRNA: long noncoding RNAs (> 200bp)
  #
  Feature_matrix$lncRNA  <- mod_gr %over% exbytx_nc[lnc_idx]
  i  = Speak(i, "lncRNA",
             "Dummy variable",
             paste0("Site overlapping with lncRNA (ncRNA > 200bp)."))

  # - HK_genes: mapped to house keeping genes, such as defined by paper below.
  #
  if (!is.null(hk_genes_id)) {
    Feature_matrix$HK_genes <-
      mod_gr %over% exbg_txdb[names(exbg_txdb) %in% hk_genes_id]
    i  = Speak(i, "HK_genes",
               "Dummy variable",
               paste0("Site overlapping with house keeping genes."))
  }

  # - miRtar_genes: mapped to microRNA targeted genes.
  #
  if (!is.null(miRtar_genes_id)) {
    Feature_matrix$miRtar_genes <-
      mod_gr %over% exbg_txdb[names(exbg_txdb) %in% miRtar_genes_id]

    i  = Speak(i, "miRtar_genes",
               "Dummy variable",
               paste0("Site overlapping with microRNA targeted genes."))
  }

  cat(paste0("Transcriptomic features generation complete.\n"))


  #Finally, mask the feature values that mapped to ambiguious genes as NA.

  if (gene_ambiguity_method == "drop_overlap") {
    Feature_matrix[mod_gr %over% removed_exbg_txdb, ] = NA
  }
  return(Feature_matrix)
}
