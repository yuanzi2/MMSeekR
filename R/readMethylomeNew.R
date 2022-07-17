#' @title Load Bis-seq data
#' @description Loading Bis-seq data from tab-delimited file or saved GRanges object
#' @import GenomicRanges
#' @importFrom dplyr left_join
#' @param fileName Text file or RDS file for GRanges object.
#' @param NNdat A dataframe of neural network score
#' @param seqLengths A named vector indicating the chromosome lengths of the genome used.
#' @param inputFormat File format. If inputFormat is set to "text" (default), the
#' argument FileName should refer to a tab-delimited text file in the
#' inputFormat: chromosome position T M, where each line stands for a CpG,
#' the position refers to the 1-base position of the methylated C (on the sense strand; if the methylated C was on antisense strand, please provide the position of G base), T is the
#' total number of reads (total counts) covering the CpG and M is the total
#' number of methylated counts. If inputFormat="GRanges", the file is assumed to be a
#' GRanges object, containing T and M as first and second data-value
#' entries, saved in rds format.
#' @return A GRanges object containing the coordinates, total (T) and methylated counts (M)
#' @export
#' @usage readMethylomeNew(fileName, NNdat, seqLengths, inputFormat = "text")
#' @examples
#' \dontrun{
#' library("MMSeekR.data")
#' library("MMSeekR")
#' data("NNscore.hg19")
#' data("hg19.seqLengths")
#' methFile <- system.file("extdata", "TCGA_BLCA_A13J_test.tab", package="MMSeekR")
#' meth <- readMethylomeNew(fileName=methFile, NNdat=NNscore.hg19, seqLengths=hg19.seqLengths)
#'}
readMethylomeNew=function(fileName, NNdat, seqLengths, inputFormat="text"){
  message("reading methylome data")
  if(inputFormat=="GRanges"){
    dat <- readRDS(fileName)
    dat=data.frame(chr=as.character(dat@seqnames), pos=as.integer(dat@ranges@start), T=dat$T, M=dat$M, stringsAsFactors = F)
  }else if(inputFormat=="text"){
    dat <- scan(fileName, what=list(chr=character(0), pos=integer(0), T=double(0), M=double(0)), sep="\t")
    dat=data.frame(dat, stringsAsFactors = F)
  }else{
    stop("unknown format")
  }
  dat=dat[dat$T>=5,]
  dat=dat[dat$chr%in%names(seqLengths),]
  dat=suppressMessages(left_join(dat, NNdat))
  meth <- GRanges(seqnames = dat$chr, ranges = IRanges(start = dat$pos, width = 1), strand = "*", T = dat$T, M = dat$M, Methylation=round(dat$M/dat$T, 6), NNscore= dat$NNScore,seqlengths = seqLengths)

  mean.cov=mean(S4Vectors::values(meth)[, 1])
  if(mean.cov < 10)
    warning(sprintf("We do not recommend the use of MMSeekR\nfor methylomes with mean coverage < 10X\n(mean coverage of CpGs with at least one read: %.1f)", mean.cov))

  meth <- meth[order(as.vector(seqnames(meth)), start(meth))]
  meth
}
