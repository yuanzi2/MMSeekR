plot3DDistributionOneChr=function(model, trainngData, name1, name2, name3, chr.sel, saveFile){
  lableList=list("Alpha score"="posterior mean of alpha", "M-value"="posterior mean of M value", "NN score" ="posterior mean of PCC")
  xStart=list("Alpha score"=0, "M-value"=(-6), "NN score"=(-1))
  xEnd=list("Alpha score"=3, "M-value"=(6), "NN score"=(1))
  pdf(saveFile, width=12, height=3.7)
  gsea.layout <- layout(matrix(1:3,nrow = 1, byrow = T))
  x=seq(xStart[[name1]], xEnd[[name1]], by=0.01)
  hist(trainngData$x[,1], probability=TRUE, breaks=30, xlab=sprintf(paste0(lableList[[name1]], " (%s)"), chr.sel), main="");
  lines(x, dnorm(x, mean=model$parms.emission$mu[[1]][[1]], sd=sqrt(model$parms.emission$sigma[[1]][[1]])), type='l', col="red");
  lines(x, dnorm(x, mean=model$parms.emission$mu[[2]][[1]], sd=sqrt(model$parms.emission$sigma[[2]][[1]])), type='l', col="green");

  x=seq(xStart[[name2]], xEnd[[name2]], by=0.01)
  hist(trainngData$x[,2], probability=TRUE, breaks=30, xlab=sprintf(paste0(lableList[[name2]], " (%s)"), chr.sel), main="");
  lines(x, dnorm(x, mean=model$parms.emission$mu[[1]][[2]], sd=sqrt(model$parms.emission$sigma[[1]][[5]])), type='l', col="red");
  lines(x, dnorm(x, mean=model$parms.emission$mu[[2]][[2]], sd=sqrt(model$parms.emission$sigma[[2]][[5]])), type='l', col="green");

  x=seq(xStart[[name3]], xEnd[[name3]], by=0.01)
  hist(trainngData$x[,3], probability=TRUE, breaks=30, xlab=sprintf(paste0(lableList[[name3]], " (%s)"), chr.sel), main="");
  lines(x, dnorm(x, mean=model$parms.emission$mu[[1]][[3]], sd=sqrt(model$parms.emission$sigma[[1]][[9]])), type='l', col="red");
  lines(x, dnorm(x, mean=model$parms.emission$mu[[2]][[3]], sd=sqrt(model$parms.emission$sigma[[2]][[9]])), type='l', col="green");
  dev.off()
}
