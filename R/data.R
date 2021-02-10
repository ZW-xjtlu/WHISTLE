#' MEA RNA 2ndary structures of mm10 exons predicted by RNAfold
#'
#' The full length transcripts are extracted from mm10;
#' the transcript annotation file was downloaded from the refSeq archive-2015-07-17-14-33-26.
#' The temperature of RNAfold prediction was set at 70 degree, and
#' the gammar for MEA structure was set at 0.1. The maximum paring distance was set at 150bp.
#'
#' For the transcripts longer than 8000bp, the predictions were conducted in windows of 2000bp with steps of 1000bp.
#'
#'
#'
#' @format A GRangesList object of length 34597:
#' \describe{
#'   Each element of the GRangesList represents a transcript, exept that the first element refers to the mitochondria chromosome.
#'   The GRanges within each of the GRangesList element represent the genomic locations of the MEA RNA 2ndary structures predicted by the method above.
#' }
"Struc_mm10"

#' MEA RNA 2ndary structures of hg19 exons predicted by RNAfold
#'
#' The full length transcripts are extracted from hg19;
#' the transcript annotation file was downloaded from the refSeq archive-2015-07-17-14-32-32.
#' The temperature of RNAfold prediction was set at 70 degree, and
#' the gammar for MEA structure was set at 0.1. The maximum paring distance was set at 150bp.
#'
#' For the transcripts longer than 8000bp, the predictions were conducted in windows of 2000bp with steps of 1000bp.
#'
#' @format A GRangesList object of length 53350:
#' \describe{
#'   Each element of the GRangesList represents a transcript, exept that the first element refers to the mitochondria chromosome.
#'   The GRanges within each of the GRangesList element represent the genomic locations of the MEA RNA 2ndary structures predicted by the method above.
#' }
"Struc_hg19"


