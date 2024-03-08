
### to obtain preliminary (single within one or two way) estimate of the system equations

preliminaryEstimate<-function(modelFrame=modelFrame, method=method){

  data=modelFrame$sysdata
  f1w<-matrix(NA,nrow=dim(data)[1],ncol=modelFrame$neq)
  f2w<-matrix(NA,nrow=dim(data)[1],ncol=modelFrame$neq)
  mi_f1w<-matrix(NA,nrow=modelFrame$nind,ncol=modelFrame$neq)
  mi_f2w<-matrix(NA,nrow=modelFrame$nind,ncol=modelFrame$neq)
  mt_f2w<-matrix(NA,nrow=length(modelFrame$nt),ncol=modelFrame$neq)
  m_f1w<-matrix(NA,nrow=1,ncol=modelFrame$neq)
  m_f2w<-matrix(NA,nrow=1,ncol=modelFrame$neq)
  Ylist<-list()
  reglist<-list()
  reglist2<-list()


  for(i in 1:modelFrame$neq){
    if(method=="1wayWB" || method=="2wayWB" || method=="2wayQUE"){
      fe1<-plm::plm(modelFrame$eqlist[[i]], data=data,
                    model="within")

      e1w <- fe1$model[,1]-fe1$coefficients%*%t(as.matrix(fe1$model[,-1]))
      f1w[,i]<-e1w-mean(e1w)
      mi_f1w[,i] <-stats::aggregate(x = f1w[,i],
                             by = list(plm::index(data)[,1]),
                             FUN = mean)$x
      m_f1w[,i] <- mean(f1w[,i])
      Ylist[[i]]<-as.matrix(fe1$model)[,1]
      reglist[[i]]<-as.matrix(fe1$model)[,-1]
      reglist2[[i]]<-fe1$model[-1]
    }

    if(method=="2wayWB" || method=="2wayQUE"){
      fe2<- plm::plm(modelFrame$eqlist[[i]],
                     data=data,
                     model="within",
                     effect="twoways")

      e2w <- fe2$model[,1]-fe2$coefficients%*%t(as.matrix(fe2$model[,-1]))
      f2w[,i]<-e2w-mean(e2w)
      mi_f2w[,i] <-stats::aggregate(x = f2w[,i],
                             by = list(plm::index(data)[,1]),
                             FUN = mean)$x

      mt_f2w[,i] <-stats::aggregate(x = f2w[,i],
                             by = list(plm::index(data)[,2]),
                             FUN = mean)$x
      m_f2w[,i] <- mean(f2w[,i])
    }}


  regnames<-NA
  for(m in 1:modelFrame$neq){
    A<-names(reglist2[[m]])
    regnames<-c(regnames,"const", A)
  }

  regnames <- regnames[-1]
  final_regnames<-regnames

  if(is.null(modelFrame$constr)){
    final_regnames=final_regnames
  }else{
    if(length(dim(modelFrame$constr)) == 0){
      z=modelFrame$constr[2]
      final_regnames<-final_regnames[-z]
    }else{
    for(i in dim(modelFrame$constr)[1]:1){
      z<-modelFrame$constr[i,2]
      final_regnames<-final_regnames[-z]
    }}}


  return(list("f1w"=f1w,
              "f2w"=f2w,
              "mi_f1w"=mi_f1w,
              "mi_f2w"=mi_f2w,
              "mt_f2w"=mt_f2w,
              "m_f1w"=m_f1w,
              "m_f2w"=m_f2w,
              "reglist"=reglist,
              "reglist2"=reglist2,
              "Ylist"=Ylist,
              "regnames"=regnames,
              "final_regnames"=final_regnames))
}
