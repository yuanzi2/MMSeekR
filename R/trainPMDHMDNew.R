#' @title Training PMD and HMD models
#' @description Train a two-state Hidden Markov Model (HMM) with Gaussian emissions. Some code was adapted from the MethylSeekR package (https://git.bioconductor.org/packages/MethylSeekR).
#' @import MethylSeekR zoocat
#' @param m A GRanges object containing the coordinates, total (T) and methylated counts (M)
#' @param chr.sel Chromosome on which HMM should be trained. Must be one of the sequence levels of m. (default chr16)
#' @param nCGbin The number of CpGs in each sliding window used to calculate PCC correlation and Mvalue (default 201). The default is highly recommended.
#' @param num.cores The number of cores used for the calculations (default 1)
#' @param pdfFilename Name of the pdf file in which the figure is saved.
#' @return A list of the parameters of the HMM model
#' @export
#' @usage trainPMDHMDNew(m, chr.sel, nCGbin=201, num.cores=1, pdfFilename)
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
#' }
trainPMDHMDNew <-function(m, chr.sel, nCGbin=201, num.cores=1, pdfFilename){
  ##calculate alpha values
  message("training PMD-HMM on chromosome ", chr.sel)
  indx <- as.character(seqnames(m))==chr.sel
  if(sum(indx)<nCGbin)
    stop(sprintf("Error: less than %d covered CpGs on chromosome %s", nCGbin, chr.sel))
  methTemp=as.data.frame(m[indx])
  T=as.numeric(methTemp$T)
  M=as.numeric(methTemp$M)
  alphaScore <- calculateAlphaDistr(M, T, nCGbin, num.cores)
  methTemp$alphaScore=alphaScore

  ##calculte the Mvalue and PCC correlation
  methTemp2=methTemp[!is.na(methTemp$NNscore),]
  NNScore <- as.numeric(methTemp2$NNscore)
  methylation=as.numeric(methTemp2$Methylation)
  methylationMean=as.vector(runmean(Rle(methylation), k = nCGbin, na.rm = TRUE, endrule = "constant"))
  methTemp2$methylationMean=methylationMean
  methTemp2$MValue=log2((methTemp2$methylationMean+0.01)/(1-methTemp2$methylationMean+0.01))
  corResult=as.vector(rollcor(NNScore, methylation, width = nCGbin,show = F, use="na.or.complete"))
  corResult=c(rep(corResult[1],(nCGbin-1)/2), corResult, rep(corResult[length(corResult)],(nCGbin-1)/2))
  methTemp2$cor=corResult

  ##calculate the model cutoff
  methTemp3=methTemp2[order(methTemp2$Methylation),]
  bottom10Per=mean(methTemp3[1:round(nrow(methTemp3)/10,0),]$Methylation)
  # top10Per=mean(methTemp3[(nrow(methTemp3)-(round(nrow(methTemp3)/10,0))+1):nrow(methTemp3),]$Methylation)
  rm(methTemp3)

  ##plot the gaussian curve
  hmmfit2=function(x,start.val,mstep){
    out<-tryCatch(
      {
        hmmfit(x, start.val, mstep)
      },
      error=function(cond){
        message(cond)
        message("")
        return(NA)
      }
    )
    return(out)
  }
  chooseModel=""
  #if(bottom10Per<0.027&top10Per>0.957){
  if(bottom10Per<0.025){
    message("Chose Model 3D")
    trainAlphaMValueNNscore <-list(x=methTemp2[,colnames(methTemp2)%in%c("alphaScore", "MValue","cor")], N=nrow(methTemp2))
    startval = hmmfit2(trainAlphaMValueNNscore, initial.model$model3D, mstep=mstep.mvnorm)
    if(sum(is.na(startval))>0){
      stop(paste0(chr.sel, " failed to train the prediction model. Maybe try other chomosome for trainning"))
    }else{
      startval=startval$model
    }
    chooseModel="Model3D"
    plot3DDistributionOneChr(startval, trainAlphaMValueNNscore, "Alpha score", "M-value", "NN score", chr.sel,pdfFilename)
  }else{
    message("Chose Model 2D")
    trainAlphaNNscore <-list(x=methTemp2[,colnames(methTemp2)%in%c("alphaScore", "cor")], N=nrow(methTemp2))
    startval = hmmfit2(trainAlphaNNscore, initial.model$model2D, mstep=mstep.mvnorm)
    if(sum(is.na(startval))>0){
      stop(paste0(chr.sel, " failed to train the prediction model. Maybe try other chomosome for trainning"))
    }else{
      startval=startval$model
    }
    chooseModel="Model2D"
    plot2DDistributionOneChr(startval, trainAlphaNNscore,  "Alpha score", "NN score", chr.sel, pdfFilename)
  }
  startval$chooseModel=chooseModel
  startval
}
