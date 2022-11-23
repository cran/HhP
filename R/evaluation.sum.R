#' Hierarchical Heterogeneity Regression Analysis.
#'
#' @usage evaluation.sum(n,q,p,admmres, abic.n, admmres2, Beta0, bic.var)
#'
#' @description The main function for Transfer learning for tensor graphical models.
#' @param n The sample size.
#' @param q The dimension of type 1 features.
#' @param p The dimension of type 2 features.
#' @param admmres The results corresponding to lambda1.
#' @param abic.n The BIC values.
#' @param admmres2 The results corresponding to lambda1.
#' @param Beta0 The true values of beta.
#' @param bic.var The BIC values.
#'
#' @return A result list including: evaluating indicator
#'
#' @export


evaluation.sum  =  function(n,q,p,admmres, abic.n, admmres2, Beta0, bic.var)
{


  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: evaluation.tuning
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##    Evaluating performances of proposed methods under all tuning parameters,
  ##    and selecting the optional results.
  ## -----------------------------------------------------------------------------------------------------------------

  sample.index  =  utils::combn(n,2)
  a.matrix.main  =  kronecker(Matrix(t(apply(sample.index,2,dMatrixFun,n=n)),sparse = T),bdiag(diag(q)))
  a.matrix.GE  =  kronecker(Matrix(t(apply(sample.index,2,dMatrixFun,n=n)),sparse = T),bdiag(diag(p)))

  L1  =  length(admmres)
  L2  =  length(admmres2)
  err.all  =  matrix(0,L1,6)
  for (l1 in 1:L1) {
    ad  =  admmres[[l1]]
    reseva  =  evaluation(n,q, p, ad$beta.gfl,ad$eta.gfl,Beta0,
                         sample.index, a.matrix.main, a.matrix.GE)
    RI.main  =  reseva$RI.main
    RI.GE  =  reseva$RI.GE
    err.all[l1,]  =  c(RI.main$gr.num,1-RI.main$group.error[1],reseva$error.main,RI.GE$gr.num,1-RI.GE$group.error[1],reseva$error.GE)
  }

  err.all2  =  matrix(0,L2,6)
  for (l2 in 1:L2) {
    ad  =  admmres2[[l2]]
    reseva  =  evaluation(n,q, p, ad$beta.gfl,ad$eta.gfl,Beta0,
                         sample.index, a.matrix.main, a.matrix.GE)
    RI.main  =  reseva$RI.main
    RI.GE  =  reseva$RI.GE
    err.all2[l2,]  =  c(RI.main$gr.num,1-RI.main$group.error[1],reseva$error.main,RI.GE$gr.num,1-RI.GE$group.error[1],reseva$error.GE)
  }

  gr.main.n  =  err.all[abic.n,1]
  n.big  =  which(err.all[,4] > gr.main.n)
  if(length(n.big) > 0){
    bic.part  =  bic.var[n.big]
    err.part  =  err.all[n.big,]
    if(length(n.big) == 1){err.part  =  t(as.matrix(err.all[n.big,]))}
    abic.part  =  which(bic.part == min(bic.part[!is.na(bic.part)]))[1]
    err.s.big  =  err.part[abic.part,]
  }
  if(length(n.big) == 0){
    err.s.big  =  err.all[abic.n,]
    err.s.big[4:6]  =  0
  }

  err.s.big = as.data.frame(t(err.s.big))
  names(err.s.big)  =  c(paste0("subgroup-",c("K","RI","MSE")),paste0("sub-subgroup-",c("K","RI","MSE")))

  result = list(err.s=err.s.big)

  return(result)

}
