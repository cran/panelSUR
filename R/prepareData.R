
### to prepare data that have to be used

prepareData<-function(data=data, restrictions=NULL, eqlist=eqlist)
{
  nconstr<-length(restrictions)
  neq<-length(eqlist)

  varlist<-list()

  for(i in 1:neq){
    varlist[[i]]<-all.vars(eqlist[[i]])
    nvar<-length(unique(unlist(varlist,use.names=FALSE)))-neq
  }

  ncoeff<-sum(sapply(varlist,length))-nconstr

  sumreg<-list()
  sumreg0<-1
  sumreg[[1]]<-sumreg0

  for(m in 2:neq){
    sumreg1<-length(all.vars(eqlist[[m-1]]))
    sumreg2<-sumreg0+sumreg1
    sumreg[[m]]<-sumreg2
    sumreg0<-sumreg2
  }

  sysdata<-subset(data, select=as.vector(unique(unlist(varlist)))) #selection of the system variables
  sysdata<-sysdata[stats::complete.cases(sysdata),]

  ####

  psur<-table(plm::index(sysdata)[,1])                #number of times each individual is observed

  sysdata$PSUR <- NA
  for(level in levels(plm::index(sysdata)[,1])) {
    subset_df <- subset(sysdata, plm::index(sysdata)[,1]==level)
    for(j in 1:nrow(subset_df)){
      subset_df$PSUR[j]<-psur[which(names(psur)==level)]}
    sysdata[plm::index(sysdata)[,1]==level, ] <- subset_df
  }

  ####


  nind<-length(unique(plm::index(sysdata)[,1]))       #number of individuals
  nt<-table(plm::index(sysdata)[,2])                  #number of individuals per group
  tmax<-length(nt)


  sumTi<-0
  for(i in 1:nind){
    Ti<-psur[i]^2
    sumTiold=sumTi
    sumTi=sumTiold+Ti
  }

  sumNt<-0
  for(t in 1:tmax){
    Nt<-nt[t]^2
    sumNtold=sumNt
    sumNt=sumNtold+Nt
  }

  vectorTi<-as.vector(table(plm::index(sysdata)[,1]))

  infoSample <- plm::pdim(sysdata)

  psurmax <- max(infoSample$Tint$Ti)



### to convert the constraint/s into the program notation

  if(is.null(restrictions)){
    constr=NULL
  }else{
    constr<-matrix(NA, ncol=2, nrow=nconstr)
    restrict<-strsplit(restrictions, "\\=")
    for(i in 1:nconstr){
      for(j in 1:2){
        restr<-strsplit(restrict[[i]][j],"\\$")
        eq_num <- which(sapply(eqlist, function(eq) eq == get(restr[[1]][1])))
        if(restr[[1]][2]=="const"){
          var_num=1}
        else{
          var_num<-which(formula.tools::get.vars(get(restr[[1]][1]))==restr[[1]][2])
          }
        constr[i,j]<-(sumreg[[eq_num]]-1)+var_num
      }
      constr[i,]<-sort(constr[i,])
    }
    constr<-constr[order(constr[,2]),]
    }


  return(list("eqlist"=eqlist,
              "nconstr"=nconstr,
              "neq"=neq,
              "varlist"=varlist,
              "ncoeff"=ncoeff,
              "sumreg"=sumreg,
              "nind"=nind,
              "nt"=nt,
              "psur"=psur,
              "psurmax"=psurmax,
              "tmax"=tmax,
              "sumTi"=sumTi,
              "sumNt"=sumNt,
              "vectorTi"=vectorTi,
              "infoSample"=infoSample,
              "sysdata"=sysdata,
              "constr"=constr))
}
