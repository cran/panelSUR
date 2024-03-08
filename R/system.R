
### to build and solve the system obtaining beta coefficient estimate

system<-function(modelFrame=modelFrame, firstEstimate=firstEstimate, matrices=matrices){

  data=modelFrame$sysdata
  XSXQ<-matrix(0,nrow=modelFrame$ncoeff,ncol=modelFrame$ncoeff)
  XSYQ<-matrix(0,nrow=modelFrame$ncoeff,ncol=1)

  pgroup <- length(names(table(data$PSUR)))

  for(g in 1:pgroup){
    t=as.integer(names(table(data$PSUR))[g])
    pM=t*modelFrame$neq
    Ip<-diag(t)
    Jp<-matrix(1/t,nrow=t,ncol=t)
    Ep<-Ip-Jp

    OpInvQ<-fastmatrix::kronecker.prod(Ep,solve(matrices$Sigma_u+matrices$Sigma_nu))+fastmatrix::kronecker.prod(Jp,solve(matrices$Sigma_u+t*matrices$Sigma_mu+matrices$Sigma_nu))

    group <- names(modelFrame$psur[which(modelFrame$psur==t)])
    nlevel<-length(group)

    for(i in 1:nlevel){
      Ysur<-matrix(0,nrow=pM,ncol=1)
      Xsur<-matrix(0,nrow=pM,ncol=modelFrame$ncoeff+modelFrame$nconstr) #ncol= # parameters excluding intercepts
      for(j in 1:t){
        for(m in 1:modelFrame$neq){
          k<-names(data$PSUR[which(plm::index(data)[,1]==group[i])])[j]
          Ysur[(j-1)*modelFrame$neq+m,1]<-firstEstimate$Ylist[[m]][which(names(firstEstimate$Ylist[[m]])==k)]
          Xsur[(j-1)*modelFrame$neq+m,modelFrame$sumreg[[m]]]<-1
          for(s in 1:length(firstEstimate$reglist2[[m]])){
            Xsur[(j-1)*modelFrame$neq+m,modelFrame$sumreg[[m]]+s]<-firstEstimate$reglist2[[m]][,s][which(names(firstEstimate$Ylist[[m]])==k)]
          }
        }
      }


      #to incorporate the restrictions into the regressors' matrix
      if(is.null(modelFrame$constr)){
        Xsur=Xsur
      }else{
        if(length(dim(modelFrame$constr)) == 0){
          z=modelFrame$constr[1]
          b=modelFrame$constr[2]
          Xsur[,z]<-Xsur[,z]+Xsur[,b]
          Xsur<-Xsur[,-b]
          }else{
        for(r in 1:dim(modelFrame$constr)[1]){
          z<-modelFrame$constr[r,1]
          b<-modelFrame$constr[r,2]
          Xsur[,z]<-Xsur[,z]+Xsur[,b]
        }
        for(i in dim(modelFrame$constr)[1]:1){
          z<-modelFrame$constr[i,2]
          Xsur<-Xsur[,-z]
        }}}

      XXSQ<-t(Xsur)%*%OpInvQ%*%Xsur
      XYSQ<-t(Xsur)%*%OpInvQ%*%Ysur

      XSXQold=XSXQ
      XSXQ=XSXQold+XXSQ

      XSYQold=XSYQ
      XSYQ=XSYQold+XYSQ

    }}

  BsurQ<-solve(XSXQ)%*%XSYQ

  std_error<-sqrt(diag(solve(XSXQ)))
  t_stat<-BsurQ/std_error
  p_value<-round(stats::pt(abs(t_stat),df=dim(data)[1]-length(t_stat),lower.tail = F),5)*2


  return(list("BsurQ"=BsurQ,
              "std_error"=std_error,
              "t_stat"=t_stat,
              "p_value"=p_value))
}
