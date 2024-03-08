
### to obtain the final output table

printSUR <- function(object){

  text_tbl <- data.frame(
    Indep = object$varnames,
    Estimate = object$Estimate,
    StdErr = object$std_error,
    tstat = object$tstat,
    pvalue = object$pvalue
  )

  if(object$method=="1wayWB"){
    options(caption = "Method: One-way WB")
  }
  if(object$method=="2wayWB"){
    options(caption = "Method: Two-way WB")
  }
  if(object$method=="2wayQUE"){
    options(caption = "Method: Two-way QUE")
  }

  infoSample <- object$infoSample
  names_list <- sapply(1:object$neq, function(i) paste("R", i, " = ", round(object$Rsquared[[i]],6), sep=""))
  Rvalues <- paste(unlist(names_list), collapse = ", ")

  ### to print the final edited table

    table <- function(text_tbl){
    cat("\n")
    cat(sprintf("%-30s\n", "SUR estimation results"))
    cat(sprintf("%-30s\n", options("caption")))
    cat("\n")
    #cat(sprintf("%-50s\n", options(data)))
    print(infoSample)
    cat("\n")
    cat("=============  ==============  ============  ===========  ========\n")
    cat(" Coefficient      Estimate       Std.Error     t-value    p-value\n")
    cat("=============  ==============  ============  ===========  ========\n")
    for (i in 1:nrow(text_tbl)) {
      cat(sprintf("%-15s %12.5f  %12.5f  %11.5f  %8.5f\n",
                  text_tbl$Indep[i], text_tbl$Estimate[i], text_tbl$StdErr[i], text_tbl$tstat[i], text_tbl$pvalue[i]))
    }
    cat("=============  ==============  ============  ===========  ========\n")
    cat("\n")
    cat("Multiple R-Squared for single equation\n")
    cat(Rvalues)
    cat("\n")
    cat("\n")
  }

  table(text_tbl)

}

