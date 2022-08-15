#' @title Predict PMD and HMD segments
#' @description Predict PMD and HMD segments using the Viterbi algorithm. Some code was adapted from the MethylSeekR package (https://git.bioconductor.org/packages/MethylSeekR).
#' @import MethylSeekR zoocat GenomicRanges
#' @param m A GRanges object containing the coordinates, total (T) and methylated counts (M)
#' @param hmm.modelList A list of the parameters of the HMM model
#' @param nCGbin The number of CpGs in each sliding window used to calculate PCC correlation and Mvalue (default 201). The default is highly recommended.
#' @param num.cores The number of cores used for the calculations (default 1)
#' @return A list of the predicted results.
#' @export
#' @usage PMDviterbiSegNew(m, hmm.modelList, nCGbin=201, num.cores=1)
#' @examples
#' \dontrun{
#' library("MMSeekR.data")
#' library("MMSeekR")
#' data("NNscore.hg19")
#' data("hg19.seqLengths")
#' methFile <- system.file("extdata", "TCGA_BLCA_A13J_test.tab", package="MMSeekR")
#' if(!dir.exists("TestResult/")){dir.create("TestResult/")}
#' meth <- readMethylomeNew(fileName=methFile, NNdat=NNscore.hg19, seqLengths=hg19.seqLengths)
#' hmm.modelList=trainPMDHMDNew(meth, "chr16", 201, 1, "TestResult/TCGA_BLCA_A13J_test.multiModel.pdf")
#' y.list=PMDviterbiSegNew(meth, hmm.modelList, 201, 1)
#' }

PMDviterbiSegNew <-function(m, hmm.modelList, nCGbin=201, num.cores=1){
  message("performing viterbi segmentation")
  nCGsPerChr=table(as.character(seqnames(m)))
  chrs=names(nCGsPerChr)[nCGsPerChr>=nCGbin]
  y.list=mclapply(chrs, function(chr.sel){
    indx=as.character(seqnames(m))==chr.sel;
    methTemp=as.data.frame(m[indx])
    T=as.numeric(methTemp$T)
    M=as.numeric(methTemp$M)
    alphaScore <- calculateAlphaDistr(M, T, nCGbin, num.cores)
    methTemp$alphaScore=alphaScore
    methTemp2=methTemp[!is.na(methTemp$NNscore),]
    NNScore <- as.numeric(methTemp2$NNscore)
    methylation=as.numeric(methTemp2$Methylation)
    methylationMean=as.vector(runmean(Rle(methylation), k = nCGbin, na.rm = TRUE, endrule = "constant"))
    methTemp2$methylationMean=methylationMean
    methTemp2$MValue=log2((methTemp2$methylationMean+0.01)/(1-methTemp2$methylationMean+0.01))
    corResult=as.vector(rollcor(NNScore, methylation, width = nCGbin,show = F, use="na.or.complete"))
    corResult=c(rep(corResult[1],(nCGbin-1)/2), corResult, rep(corResult[length(corResult)],(nCGbin-1)/2))
    methTemp2$cor=corResult

    if(hmm.modelList$chooseModel=="Model3D"){
      trainData <-list(x=methTemp2[,colnames(methTemp2)%in%c("alphaScore", "MValue","cor")], N=nrow(methTemp2))
    }else if(hmm.modelList$chooseModel=="Model2D"){
      trainData <-list(x=methTemp2[,colnames(methTemp2)%in%c("alphaScore", "cor")], N=nrow(methTemp2))
    }
    yData=predict(hmm.modelList, trainData)
    #remove regions that are too short
    ttt=Rle(yData$s)
    min.len=101
    # first take regions that are PMD, but too short and make them nonPMD
    indx=runLength(ttt)<=min.len & runValue(ttt)==2;
    runValue(ttt)[indx]=1;
    # now vice versa
    indx=runLength(ttt)<=min.len & runValue(ttt)==1;
    runValue(ttt)[indx]=2;
    yData$s=as.vector(ttt)
    return(yData)
  }, mc.cores=num.cores);
  names(y.list)=chrs
  y.list
}
