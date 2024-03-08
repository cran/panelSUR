
### post-estimation indicators

postEstimation<-function(modelFrame=modelFrame, firstEstimate=firstEstimate, system=system){

  data=modelFrame$sysdata

  #to obtain the ordered and complete vector of coefficients

  coeff <- sapply(1:length(firstEstimate$regnames), function(i) paste(i))
  coeff_final <- coeff

  if(is.null(modelFrame$constr)){
    coeff_final=coeff
    }else{
      if(length(dim(modelFrame$constr)) == 0){
        z=modelFrame$constr[1]
        b=modelFrame$constr[2]
        coeff_final[b]<-coeff_final[z]
      }else{
    for(r in 1:dim(modelFrame$constr)[1]){
      z<-modelFrame$constr[r,1]
      b<-modelFrame$constr[r,2]
      coeff_final[b]<-coeff_final[z]
    }}
    coeff_final <- coeff_final[!duplicated(coeff_final)]}

  estRes <- cbind(coeff_final, system$BsurQ)
  colnames(estRes)<-c("coeff","value")

  beta <- merge(as.data.frame(coeff), estRes, by="coeff", all=T)
  beta=beta[order(as.integer(beta[,1])),]

  if(is.null(modelFrame$constr)){
    beta=beta
  }else{
    if(length(dim(modelFrame$constr)) == 0){
      z=modelFrame$constr[1]
      b=modelFrame$constr[2]
      beta[b,2]<-beta[z,2]
    }else{
      for(r in 1:dim(modelFrame$constr)[1]){
      z<-modelFrame$constr[r,1]
      b<-modelFrame$constr[r,2]
      beta[b,2]<-beta[z,2]
    }}}

  beta <- as.matrix(as.numeric(beta[,2]))

  #to obtain the Rsquared values for each single equation

  Rsquared <- list()

  for(i in 1:modelFrame$neq){
    u <- matrix(1,nrow=dim(data)[1], ncol=1)
    newX <- cbind(u,firstEstimate$reglist[[i]])
    t<-length(modelFrame$varlist[[i]])
    Yhat<-newX%*%beta[modelFrame$sumreg[[i]]:(modelFrame$sumreg[[i]]+t-1)]
    error <- firstEstimate$Ylist[[i]]-Yhat
    Rsquared[[i]] <- (t(error)%*%error)/(t(firstEstimate$Ylist[[i]]-mean(firstEstimate$Ylist[[i]]))%*%(firstEstimate$Ylist[[i]]-mean(firstEstimate$Ylist[[i]])))
  }

  return(list("Rsquared"=Rsquared))
}
