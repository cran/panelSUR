

### final command

SURest<-function(data=data,eqlist=eqlist, restrictions=NULL, method="1wayWB"){

  ## checking the correct specification of equations
  for(i in 1:length(eqlist)){
    if (!is.null(attr(terms(eqlist[[i]]), "intercept"))) {
      if (!attr(terms(eqlist[[i]]), "intercept")) {
        stop("The constant term must be included in each equation")
      }
    } else {
      stop("The constant term must be included in each equation.")
    }}


  ## checking eqlist
  if(!inherits(eqlist, "list" ) ){
    stop( "eqlist must be of the class 'list'" )
  }

  ## checking eqlist contents
  for(i in 1:length(eqlist)){
    if(!inherits(eqlist[[i]], "formula" )){
      stop( "eqlist elements must be of the class 'formula'" )
    }}

  ## checking data type
  if(!(class( data )[1] == "pdata.frame" ) ){
    stop( "Data must be of the class 'pdata.frame'" )
  }

  ## checking argument 'method'
  if(!( method %in% c("1wayWB", "2wayWB", "2wayQUE") ) ){
    stop( "The method must be '1wayWB', '2wayWB', '2wayQUE'" )
  }


  modelFrame<-prepareData(data=data, restrictions=restrictions, eqlist=eqlist)

  firstEstimate<-preliminaryEstimate(modelFrame=modelFrame, method=method)

  matrices<-obtainSigmas(modelFrame=modelFrame, firstEstimate=firstEstimate, method=method)

  system<-system(modelFrame=modelFrame, firstEstimate=firstEstimate, matrices=matrices)

  postEst <- postEstimation(modelFrame=modelFrame, firstEstimate=firstEstimate, system=system)

  return(list("Sigma_mu"=matrices$Sigma_mu,
              "Sigma_nu"=matrices$Sigma_nu,
              "Sigma_u"=matrices$Sigma_u,
              "varnames"=firstEstimate$final_regnames,
              "Estimate"=system$BsurQ,
              "std_error"=system$std_error,
              "tstat" = system$t_stat,
              "pvalue"=system$p_value,
              "infoSample"=modelFrame$infoSample,
              "neq"=modelFrame$neq,
              "Rsquared"=postEst$Rsquared,
              "method"=method))
}
