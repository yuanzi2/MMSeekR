#' @title Predict PMDs using multi-model(MMSeekR) or singel-model (MethylSeekR origin model)
#' @description PMD prediction using  multi-model(MMSeekR) or singel-model (MethylSeekR origin model)
#' @import GenomicRanges MethylSeekR zoocat
#' @importFrom dplyr left_join
#' @importFrom valr bed_subtract
#' @param genomeVersion genome version. hg19 or hg38 for human.
#' @param fileName Methylome file. Text file or RDS file for GRanges object.
#' @param outputIndex Ouput file index, which contains outputPath and outputPrefix, like "outputPath/outputPrefix"
#' @param trainChr Chromosome on which HMM should be trained. Must be one of the sequence levels in Methylome file. (default chr16).
#' If you want to run singel-model (MethylSeekR origin model) using runMethylSeekR=T, you can set trainChr as "chr16" (MethylSeekR default parameters) or set useMethylSeekRDefaultParameters=TRUE
#' @param inputFormat File format. If format is set to "text" (default), the
#' argument FileName should refer to a tab-delimited text file in the
#' format: chromosome position T M, where each line stands for a CpG,
#' the position refers to the 1-base position of the methylated C (on the sense strand; if the methylated C was on antisense strand, please provide the position of G base), T is the
#' total number of reads (total counts) covering the CpG and M is the total
#' number of methylated counts. If format="GRanges", the file is assumed to be a
#' GRanges object, containing T and M as first and second data-value
#' entries, saved in rds format.
#' @param num.cores The number of cores used for the calculations (default 1)
#' @param nCGbin The number of CpGs in each sliding window used to calculate PCC correlation and Mvalue. This number should be an odd number. The default (201) is highly recommended for multi-model.
#' If you want to run singel-model (MethylSeekR origin model) using runMethylSeekR=T,  you can set nCGbin as "101" (MethylSeekR default parameters) or set useMethylSeekRDefaultParameters=TRUE
#' @param rmSNP Whether to remove SNPs (default FALSE). If yes, SNPFile should not be NULL.
#' @param SNPFile A text file for SNP file. This file should refer to a tab-delimited text file in the format: chromosome position, where each line stands for a SNP. It only works when rmSNP=T.
#' @param runMethylSeekR Choose the singel-model (MethylSeekR origin model) to call PMDs (default FALSE).
#' @param useMethylSeekRDefaultParameters Use singel-model (MethylSeekR origin model) default parameters (tainChr="chr16" and nCGbin=101). This parameter works only when runMethylSeekR=T.
#' If you want to use singel-model (MethylSeekR origin model) to call PMDs, runMethylSeekR=T is recommended.
#' @return No return value
#' @export
#' @usage runPMDs(genomeVersion, fileName, outputIndex, trainChr="chr16", inputFormat="text",
#'                num.cores = 1, nCGbin = 201, rmSNP=FALSE, SNPFile=NULL, runMethylSeekR=FALSE,
#'                useMethylSeekRDefaultParameters=FALSE)
#' @examples
#' library("MMSeekR.data")
#' library("MMSeekR")
#' methFile <- system.file("extdata", "TCGA_BLCA_A13J_test.tab", package="MMSeekR")
#' ## Predict PMDs using multi-model(new model) without removing SNPs
#' runPMDs("hg19", methFile, "TestResult/TCGA_BLCA_A13J_test")
#' @examples
#' \dontrun{
#' ## Predict PMDs using multi-model(MMSeekR) with removing SNPs
#' SNPFile <- system.file("extdata", "common_chr16_22_20180418.tab", package="MMSeekR")
#' runPMDs("hg19", methFile, "TestResult/TCGA_BLCA_A13J_test", rmSNP=TRUE, SNPFile=SNPFile)
#' }
#'
#' @examples
#' ## Predict PMDs using singel-model (MethylSeekR origin model) with MethylSeekR default parameters and without removing SNPs
#' runPMDs("hg19", methFile, "TestResult/TCGA_BLCA_A13J_testDefault",
#'          runMethylSeekR=TRUE, useMethylSeekRDefaultParameters=TRUE)
#' @examples
#' \dontrun{
#' ## Predict PMDs using singel-model (MethylSeekR origin model) with MethylSeekR default parameters and removing SNPs
#' SNPFile <- system.file("extdata", "common_chr16_22_20180418.tab", package="MMSeekR")
#' runPMDs("hg19", methFile, "TestResult/TCGA_BLCA_A13J_testDefault", runMethylSeekR=TRUE,
#'          useMethylSeekRDefaultParameters=T, rmSNP=T, SNPFile=SNPFile)
#' }
#'
runPMDs=function(genomeVersion, fileName, outputIndex, trainChr="chr16", inputFormat="text", num.cores = 1, nCGbin = 201, rmSNP=FALSE, SNPFile=NULL, runMethylSeekR=FALSE, useMethylSeekRDefaultParameters=FALSE){
  if(runMethylSeekR){
    if(useMethylSeekRDefaultParameters){
      message("Runing single-model (MethylSeekR origin model)")
      message("Using the default parameters of MethylSeekR: taining chromosome=chr22, nCGbin=101")
      runSingleModel(genomeVersion, fileName, outputIndex, trainChr="chr22", inputFormat=inputFormat, num.cores=num.cores, nCGbin=101, rmSNP=rmSNP, SNPFile=SNPFile)
    }else{
      message("Runing single-model (MethylSeekR origin model)")
      message(paste0("Using the user's parameters: taining chromosome=", trainChr, ", nCGbin=", nCGbin))
      runSingleModel(genomeVersion, fileName, outputIndex, trainChr=trainChr, inputFormat=inputFormat, num.cores=num.cores, nCGbin=nCGbin, rmSNP=rmSNP, SNPFile=SNPFile)
    }
  }else{
    message("Runing multi-model")
    message(paste0("Using the user's parameters: taining chromosome=", trainChr, ", nCGbin=", nCGbin))
    runMultiModel(genomeVersion, fileName, outputIndex, trainChr=trainChr, inputFormat=inputFormat, num.cores = num.cores, nCGbin = nCGbin, rmSNP=rmSNP, SNPFile=SNPFile)
  }
}
