#' @title Remove blackList and save PMD segments
#' @description Remove blackList and save PMD segments in rds format and a tab-delimited file. Some code was adapted from the MethylSeekR package (https://git.bioconductor.org/packages/MethylSeekR).
#' @import GenomicRanges tibble
#' @importFrom valr bed_subtract
#' @param seg GRanges object.
#' @param blackList A tibble of blackList
#' @param GRangesFilename Filename of the GRanges object.
#' @param TableFilename Filename of the PMD table.
#' @return No return value
#' @export
#' @usage savePMDSegRmBlackListNew(seg, blackList, GRangesFilename = NULL, TableFilename = NULL)
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
#' savePMDSegRmBlackListNew(segments, hg19.blackList, "TestResult/TCGA_BLCA_A13J_test.multiModel.rds", "TestResult/TCGA_BLCA_A13J_test.multiModel.PMDs.bed")
#'}
savePMDSegRmBlackListNew <- function(seg, blackList, GRangesFilename = NULL, TableFilename = NULL){
  seg=tibble(chrom=as.character(seqnames(seg)), start=as.integer(start(seg)), end=as.integer(end(seg)), type=seg$type)
  seg=bed_subtract(seg, blackList)
  seg=GRanges(seqnames = seg$chrom, IRanges(start=seg$start, end=seg$end), strand = "*", type=seg$type)
  # save as GRanges object
  if(!is.null(GRangesFilename))
    saveRDS(seg, GRangesFilename)
  # save as tab-delimited table
  if(!is.null(TableFilename)){
    indx=S4Vectors::values(seg)$type=="PMD"
    PMDsResult=data.frame(chr=as.character(seqnames(seg))[indx], start=start(seg)[indx], end=end(seg)[indx])
    write.table(PMDsResult, file=TableFilename, quote=FALSE, sep="\t", row.names=FALSE, col.names = F)
  }
}
