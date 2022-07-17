#' @title Create an GRanges Object containing the PMD segmentation.
#' @import GenomicRanges
#' @description Create an GRanges Object containing the PMD segmentation
#' @param m A GRanges object containing the coordinates, total (T) and methylated counts (M)
#' @param y.list A list of the predicted results.
#' @param num.cores The number of cores used for the calculations (default 1)
#' @param seqLengths A named vector indicating the chromosome lengths of the genome used.
#' @return A GRanges object containing the PMD segmentation.
#' @export
#' @usage createGRangesObjectPMDSegNew(m, y.list, num.cores=1, seqLengths)
#' @examples
#' \dontrun{
#' library("MMSeekR.data")
#' library("MMSeekR")
#' data("NNscore.hg19")
#' data("hg19.seqLengths")
#' data("hg19.blackList")
#' methFile <- system.file("extdata", "TCGA_BLCA_A13J_test.tab", package="MMSeekR")
#' meth <- readMethylomeNew(fileName=methFile, NNdat=NNscore.hg19, seqLengths=hg19.seqLengths)
#' if(!dir.exists("TestResult/")){dir.create("TestResult/")}
#' hmm.modelList=trainPMDHMDNew(meth, "chr16", 201, 1, "TestResult/TCGA_BLCA_A13J_test.multiModel.pdf")
#' y.list=PMDviterbiSegNew(meth, hmm.modelList, 201, 1)
#' segments = createGRangesObjectPMDSegNew(meth, y.list, 1, hg19.seqLengths)
#' }
#'
createGRangesObjectPMDSegNew <-function(m, y.list, num.cores=1, seqLengths){
  message("creating GRanges object")
  chrs=names(y.list)
  segList=mclapply(chrs, function(chr.sel) {
    indx=as.character(seqnames(m))==chr.sel
    n <- sum(indx)
    methTemp=m[indx]
    index2=!is.na(S4Vectors::values(methTemp)$NNscore)
    n <- sum(index2)
    methTemp2=methTemp[index2]

    mids <- round(0.5*(start(methTemp2)[-length(methTemp2)] + start(methTemp2)[-1]))
    segCG <- Rle(y.list[[chr.sel]]$s)
    segNt <- Rle(lengths=c(diff(c(1,mids)),seqLengths[chr.sel]-mids[n-1]+1), values=y.list[[chr.sel]]$s)
    segChr <- GRanges(seqnames=chr.sel, IRanges(start=c(1,cumsum(runLength(segNt))[-nrun(segNt)]+1),
                                                end=cumsum(runLength(segNt))), strand="*", type=c("notPMD", "PMD")[runValue(segNt)], nCG=runLength(segCG), seqlengths=seqLengths)
    segChr
  }, mc.cores=num.cores);

  segments <- do.call(c, unname(segList))
  segments
}
