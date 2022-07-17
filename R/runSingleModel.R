#' @title Predict PMDs using single-model (the origin model of MethylSeekR)
#' @description PMD prediction using single-model
#' @import MethylSeekR GenomicRanges
#' @importFrom valr bed_subtract
#' @param genomeVersion genome version. hg19 or hg38 for human.
#' @param fileName Methylome file. Text file or RDS file for GRanges object.
#' @param outputIndex Ouput file index, which contains outputPath and outputPrefix, like "outputPath/outputPrefix"
#' @param trainChr Chromosome on which HMM should be trained. Must be one of the sequence levels in Methylome file. (default chr22)
#' @param inputFormat File format. If format is set to "text" (default), the
#' argument FileName should refer to a tab-delimited text file in the
#' format: chromosome position T M, where each line stands for a CpG,
#' the position refers to the 1-base position of the methylated C (on the sense strand; if the methylated C was on antisense strand, please provide the position of G base), T is the
#' total number of reads (total counts) covering the CpG and M is the total
#' number of methylated counts. If format="GRanges", the file is assumed to be a
#' GRanges object, containing T and M as first and second data-value
#' entries, saved in rds format.
#' @param num.cores The number of cores used for the calculations (default 1)
#' @param nCGbin The number of CpGs in each sliding window used to calculate alpha (default 101). The default is highly recommended. This number should be an odd number
#' @param rmSNP Whether to remove SNPs (default FALSE). If yes, SNPFile should not be NULL.
#' @param SNPFile A text file for SNP file. This file should refer to a tab-delimited text file in the format: chromosome position, where each line stands for a SNP. It only works when rmSNP=T.
#' @return No return value
#' @export
#' @usage runSingleModel(genomeVersion, fileName, outputIndex, trainChr="chr22", inputFormat="text", num.cores = 1, nCGbin = 101, rmSNP=FALSE, SNPFile=NULL)
#' @examples
#' \dontrun{
#' library("MMSeekR.data")
#' library("MMSeekR")
#' methFile <- system.file("extdata", "TCGA_BLCA_A13J_test.tab", package="MMSeekR.data")
#' ## Do not remove SNPs
#' runSingleModel("hg19", methFile, "TestResults/TCGA_BLCA_A13J_test")
#' ## Remove SNPs
#' SNPFile <- system.file("extdata", "common_chr16_22_20180418.tab", package="MMSeekR")
#' runSingleModel("hg19", methFile, "TestResults/TCGA_BLCA_A13J_test", rmSNP=TRUE, SNPFile=SNPFile)
#' }

runSingleModel=function(genomeVersion, fileName, outputIndex, trainChr="chr22", inputFormat="text", num.cores = 1, nCGbin = 101, rmSNP=FALSE, SNPFile=NULL){
  if(nCGbin<101){
    stop("nCGbin is recommended >= 101")
  }else if(nCGbin%%2==0){
    stop("nCGbin should be a odd number")
  }

  ###check whether the outputPath exist. If not, create the outputPath
  if(!dir.exists(dirname(outputIndex))){
    dir.create(dirname(outputIndex))
  }

  if(genomeVersion=="hg38"){
    data(list = c("hg38.blackList", "hg38.seqLengths"), package = "MMSeekR.data")
    seqLengths=hg38.seqLengths
    blackList=hg38.blackList
  }else if(genomeVersion=="hg19"){
    data(list = c("hg19.blackList", "hg19.seqLengths"), package = "MMSeekR.data")
    seqLengths=hg19.seqLengths
    blackList=hg19.blackList
  }else{
    stop("Currently only support hg19, hg38")
  }
  m=MethylSeekR::readMethylome(fileName, seqLengths, inputFormat)
  seqLengths=seqLengths[names(seqLengths)%in%as.character(unique(m@seqnames))]
  m=m[m$T>=5]
  if(rmSNP){
    if(is.null(SNPFile)){
      stop("Please provide SNP file for SNPFile parameter!")
    }else{
      snps.gr <- MethylSeekR::readSNPTable(FileName=SNPFile, seqLengths=seqLengths)
      m <- MethylSeekR::removeSNPs(m, snps.gr)
    }
  }else{
    message("Do not remove SNPs")
  }
  segs=MethylSeekR::segmentPMDs(m, trainChr, pdfFilename = paste0(outputIndex, ".singleModel.pdf"), seqLengths, num.cores = num.cores, nCGbin = nCGbin)
  savePMDSegRmBlackListNew(segs, blackList, paste0(outputIndex, ".singleModel.rds"), paste0(outputIndex, ".singleModel.PMDs.bed"))
}



