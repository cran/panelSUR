
### to obtain the errors' variance-covariance matrices

obtainSigmas<-function(modelFrame=modelFrame,firstEstimate=firstEstimate,method=method){

  data=modelFrame$sysdata
  fBc1w_1<-matrix(NA,nrow=modelFrame$neq,ncol=1)
  fBc1w_2<-matrix(NA,nrow=modelFrame$neq,ncol=modelFrame$neq)
  Bc1w<-matrix(0,nrow=modelFrame$neq,ncol=modelFrame$neq)
  fBc2w_1<-matrix(NA,nrow=modelFrame$neq,ncol=1)
  fBc2w_2<-matrix(NA,nrow=modelFrame$neq,ncol=modelFrame$neq)
  Bc2w<-matrix(0,nrow=modelFrame$neq,ncol=modelFrame$neq)
  fW1w_1<-matrix(NA,nrow=modelFrame$neq,ncol=1)
  fW2w_1<-matrix(NA,nrow=modelFrame$neq,ncol=1)
  fW1w_2<-matrix(NA,nrow=modelFrame$neq,ncol=modelFrame$neq)
  W1w<-matrix(0,nrow=modelFrame$neq,ncol=modelFrame$neq)
  W2w<-matrix(0,nrow=modelFrame$neq,ncol=modelFrame$neq)

  fBt2w_1<-matrix(NA,nrow=modelFrame$neq,ncol=1)
  fBt2w_2<-matrix(NA,nrow=modelFrame$neq,ncol=modelFrame$neq)
  Bt2w<-matrix(0,nrow=modelFrame$neq,ncol=modelFrame$neq)


  for(i in 1:modelFrame$nind){


    for(j in 1:modelFrame$neq){
      fBc1w_1[j]<-firstEstimate$mi_f1w[i,j]-firstEstimate$m_f1w[j]
    }

    fBc1w_2<-modelFrame$psur[i]*fBc1w_1%*%t(fBc1w_1)
    Bc1wold=Bc1w
    Bc1w=Bc1wold+fBc1w_2

    pmax<-modelFrame$psur[i]
    for(p in 1:pmax){
      k<-plm::index(data)[,2][which(plm::index(data)[,1]==unique(plm::index(data)[,1])[i])]
      t<-k[which(rank(k)==p)]
      for(j in 1:modelFrame$neq){
        fW1w_1[j]<-firstEstimate$f1w[which(plm::index(data)[,2]==t & plm::index(data)[,1]==unique(plm::index(data)[,1])[i]),j]-firstEstimate$mi_f1w[i,j]
        fW2w_1[j]<-firstEstimate$f2w[which(plm::index(data)[,2]==t & plm::index(data)[,1]==unique(plm::index(data)[,1])[i]),j]-firstEstimate$mi_f2w[i,j]-firstEstimate$mt_f2w[t,j]
      }

      fW1w_2<-fW1w_1%*%t(fW1w_1)
      W1wold=W1w
      W1w=W1wold+fW1w_2

      fW2w_2<-fW2w_1%*%t(fW2w_1)
      W2wold=W2w
      W2w=W2wold+fW2w_2
    }


    for(j in 1:modelFrame$neq){
      fBc2w_1[j]<-firstEstimate$mi_f2w[i,j]-firstEstimate$m_f2w[j]
    }

    fBc2w_2<-modelFrame$psur[i]*fBc2w_1%*%t(fBc2w_1)
    Bc2wold=Bc2w
    Bc2w=Bc2wold+fBc2w_2
  }


  if(method=="1wayWB"){
    Sigma_u<-W1w/(dim(data)[1]-modelFrame$nind)  #SigmaB_u_1w
    Sigma_mu<-(Bc1w-(modelFrame$nind-1)*Sigma_u)/(dim(data)[1]-modelFrame$sumTi/dim(data)[1]) #SigmaB_mi_1w
    Sigma_nu<-matrix(0,ncol=modelFrame$neq,nrow=modelFrame$neq)
  }


  if(method=="2wayWB"){

    for(t in 1:modelFrame$tmax){
      for(j in 1:modelFrame$neq){
        fBt2w_1[j]<-firstEstimate$mt_f2w[t,j]-firstEstimate$m_f1w[j]
      }
      fBt2w_2<-modelFrame$nt[t]*fBt2w_1%*%t(fBt2w_1)
      Bt2wold=Bt2w
      Bt2w=Bt2wold+fBt2w_2
    }

    #SigmaB_u
    Sigma_u<-W2w/(dim(data)[1]-modelFrame$nind-modelFrame$tmax)
    #SigmaB_mi
    Sigma_mu<-(Bc2w-(modelFrame$nind-1)*Sigma_u)/(dim(data)[1]-modelFrame$sumTi/dim(data)[1])
    #SigmaB_ni
    Sigma_nu<-(Bt2w-(modelFrame$tmax-1)*Sigma_u)/(dim(data)[1]-modelFrame$sumNt/dim(data)[1])
  }


  if(method=="2wayQUE"){

    DNf_femI<-list()

    for(j in 1:modelFrame$neq){
      DNf_femI[[j]]<-firstEstimate$mi_f2w[,j]*modelFrame$vectorTi
    }

    DT<-diag(modelFrame$nt)

    DTN<-matrix(0,nrow=length(modelFrame$nt), ncol=modelFrame$nind)
    DTNDNinv<-matrix(0,nrow=length(modelFrame$nt), ncol=modelFrame$nind)

    for(j in 1:dim(data)[1]){
      i=plm::index(data)[,1][j]
      t=plm::index(data)[,2][j]
      DTN[t,i]=1
      DTNDNinv[t,i]=1/modelFrame$vectorTi[i]
    }

    Q<-DT-(DTNDNinv%*%t(DTN))

    qn<-matrix(NA, nrow=modelFrame$neq, ncol=modelFrame$neq)

    for(i in 1:modelFrame$neq){
      for(j in 1:modelFrame$neq){
        qnx<-t(firstEstimate$f2w[,i])%*%firstEstimate$f2w[,j]-t(firstEstimate$mi_f2w[,i])%*%DNf_femI[[j]]-(t(firstEstimate$mt_f2w[,i])%*%DT-t(firstEstimate$mi_f2w[,i])%*%t(DTN))%*%MASS::ginv(Q)%*%(t(t(firstEstimate$mt_f2w[,j])%*%DT-t(firstEstimate$mi_f2w[,j])%*%t(DTN)))
        qn[j,i]<-qnx
      }}


    ######

    XnoCmI_list<-list()
    XnoCmT_list<-list()
    DNXnoCmI_list<-list()

    for(m in 1:modelFrame$neq){
      XnoCmI<-stats::aggregate(x = firstEstimate$reglist[[m]],
                        by = list(plm::index(data)[,1]),
                        FUN = mean)[-1]
      XnoCmI_list[[m]]<-XnoCmI

      XnoCmT<-stats::aggregate(x = firstEstimate$reglist[[m]],
                        by = list(plm::index(data)[,2]),
                        FUN = mean)[-1]
      XnoCmT_list[[m]]<-as.matrix(XnoCmT)

      DNXnoCmI<-matrix(NA,ncol=dim(firstEstimate$reglist[[m]])[2],nrow=modelFrame$nind)
      for(i in 1:modelFrame$nind){
        for(j in 1:dim(XnoCmI)[2]){
          DNXnoCmI[i,j]<- XnoCmI[i,j]%*%modelFrame$vectorTi[i]
          DNXnoCmI_list[[m]]<-DNXnoCmI
        }}}

    #####

    XQDX_list<-list()

    for(i in 1:modelFrame$neq){
      for(j in 1:modelFrame$neq){
        XQDX<-matrix(NA, nrow=dim(firstEstimate$reglist[[i]])[2], ncol=dim(firstEstimate$reglist[[j]])[2])
        XQDX<-t(firstEstimate$reglist[[i]])%*%firstEstimate$reglist[[j]]-t(XnoCmI_list[[i]])%*%DNXnoCmI_list[[j]]-(t(XnoCmT_list[[i]])%*%DT-t(XnoCmI_list[[i]])%*%t(DTN))%*%MASS::ginv(Q)%*%(t(t(XnoCmT_list[[j]])%*%DT-t(XnoCmI_list[[j]])%*%t(DTN)))
        XQDX_list[[(i-1)*modelFrame$neq+j]]<-XQDX
      }}

    Ksur<-matrix(NA, nrow=modelFrame$neq, ncol=modelFrame$neq)

    for(i in 1:modelFrame$neq){
      for(j in 1:modelFrame$neq){
        Ksur[i,j]<-matlib::tr(solve(XQDX_list[[(i-1)*modelFrame$neq+i]])%*%XQDX_list[[(i-1)*modelFrame$neq+j]]%*%solve(XQDX_list[[(j-1)*modelFrame$neq+j]])%*%XQDX_list[[(j-1)*modelFrame$neq+i]])
      }}

    Sigma_u<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)

    for(i in 1:modelFrame$neq){
      for(j in 1:modelFrame$neq){
        Sigma_u[i,j]=qn[i,j]/(dim(data)[1]-modelFrame$nind-(modelFrame$tmax-1)-Ksur[i,i]-Ksur[j,j]+Ksur[i,j]) ###c'era pmax e non tmax, ma pmax dipende da i
      }}

    qI<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)
    qT<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)
    kI<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)
    kT<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)
    k0<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)

    for(i in 1:modelFrame$neq){
      for(j in 1:modelFrame$neq){
        qI[j,i]<-t(as.vector(firstEstimate$mi_f2w[,i]))%*%DNf_femI[[j]]
        qT[j,i]<-t(as.vector(firstEstimate$mt_f2w[,i]))%*%DT%*%(as.vector(firstEstimate$mt_f2w[,j]))
        kI[j,i]<-matlib::tr(solve(XQDX_list[[1+(j-1)*(modelFrame$neq+1)]])%*%XQDX_list[[i+(j-1)*modelFrame$neq]]%*%solve(XQDX_list[[i+(i-1)*modelFrame$neq]])%*%t(XnoCmI_list[[i]])%*%DNXnoCmI_list[[j]])
        kT[j,i]<-matlib::tr(solve(XQDX_list[[1+(j-1)*(modelFrame$neq+1)]])%*%XQDX_list[[i+(j-1)*modelFrame$neq]]%*%solve(XQDX_list[[i+(i-1)*modelFrame$neq]])%*%t(XnoCmT_list[[i]])%*%DT%*%XnoCmT_list[[j]])
      }}


    XnoCsum_list<-list()

    for(i in 1:modelFrame$neq){
      XnoCsum<-as.vector(colSums(firstEstimate$reglist[[i]]))
      XnoCsum_list[[i]]<-XnoCsum
    }

    for(i in 1:modelFrame$neq){
      for(j in 1:modelFrame$neq){
        k0[j,i]<-1/dim(data)[1]*((XnoCsum_list[[j]])%*%solve(XQDX_list[[1+(j-1)*(modelFrame$neq+1)]])%*%XQDX_list[[i+(j-1)*modelFrame$neq]]%*%solve(XQDX_list[[i+(i-1)*modelFrame$neq]])%*%XnoCsum_list[[i]])
      }}


    lMI<-t(modelFrame$vectorTi)%*%modelFrame$vectorTi/dim(data)[1]
    lNI<-t(modelFrame$nt)%*%modelFrame$nt/dim(data)[1]

    Sigma_mu<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)
    Sigma_nu<-matrix(NA,ncol=modelFrame$neq,nrow=modelFrame$neq)

    for(i in 1:modelFrame$neq){
      for(j in 1:modelFrame$neq){
        Sigma_mu[i,j]=(((dim(data)[1]-lNI)*(qI[i,j]-(modelFrame$nind+kI[i,j]-k0[i,j]-1)*Sigma_u[i,j]))-((modelFrame$nind-lNI)*(qT[i,j]-(modelFrame$tmax+kT[i,j]-k0[i,j]-1)*Sigma_u[i,j])))/((dim(data)[1]-lMI)*(dim(data)[1]-lNI)-(modelFrame$nind-lNI)*(modelFrame$tmax-lMI))
        Sigma_nu[i,j]=(((dim(data)[1]-lMI)*(qT[i,j]-(modelFrame$tmax+kT[i,j]-k0[i,j]-1)*Sigma_u[i,j]))-((modelFrame$tmax-lMI)*(qI[i,j]-(modelFrame$nind+kI[i,j]-k0[i,j]-1)*Sigma_u[i,j])))/((dim(data)[1]-lMI)*(dim(data)[1]-lNI)-(modelFrame$nind-lNI)*(modelFrame$tmax-lMI))
      }}
  }


  return(list("Sigma_u"=Sigma_u,
              "Sigma_mu"=Sigma_mu,
              "Sigma_nu"=Sigma_nu))

}
