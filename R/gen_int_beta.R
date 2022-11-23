#' Hierarchical Heterogeneity Regression Analysis.
#'
#' @author Mingyang Ren, Qingzhao Zhang, Sanguo Zhang, Tingyan Zhong, Jian Huang, Shuangge Ma. Maintainer: Mingyang Ren <renmingyang17@mails.ucas.ac.cn>.
#' @references Mingyang Ren, Qingzhao Zhang, Sanguo Zhang, Tingyan Zhong, Jian Huang, Shuangge Ma. 2022. Hierarchical Cancer Heterogeneity Analysis Based On Histopathological Imaging Features. Biometrics, <DOI: 10.1111/biom.13544>.
#' @usage gen_int_beta(n, p, q, whole.data, subgroup=c(2,4),
#'                     ridge = FALSE, gr.init=10, lambda.min=0.0001)
#'
#' @description The main function for Transfer learning for tensor graphical models.
#' @param whole.data The input data analyzed (a list including the response and design matrix).
#' @param n The sample size.
#' @param q The dimension of type 1 features.
#' @param p The dimension of type 2 features.
#' @param subgroup When using fmrs to generate initial value, the initial value parameter of fmrs is given. Randomly divide this number of groups into several groups.
#' @param ridge The logical variable, whether or not to yield initial values using ridge regression.
#' @param gr.init The subgroup number of initial values using ridge regression.
#' @param lambda.min The tuning parameter using ridge regression, the default is 0.0001.
#'
#' @return A result list.
#' @export
#'
#' @import Matrix MASS fmrs
#' @importFrom methods slot
#' @importFrom utils combn
#'
#'
#' @examples
#' \donttest{
#' library(HhP)
#' library(Matrix)
#' library(MASS)
#' library(fmrs)
#' data(example.data.reg)
#' n = example.data.reg$n
#' q = example.data.reg$q
#' p = example.data.reg$p
#'
#' beta.init.list  =  gen_int_beta(n, p, q, example.data.reg)
#' beta.init  =  beta.init.list$beta.init
#' lambda  =  genelambda.obo()
#' result  =  HhP.reg(lambda, example.data.reg, n, q, p, beta.init)
#' index.list  =  evaluation.sum(n,q,p, result$admmres, result$abic.n,
#'                result$admmres2, example.data.reg$Beta0, result$bic.var)
#' index.list$err.s
#' }
#'
#'

gen_int_beta  =  function(n, p, q, whole.data, subgroup=c(2,4),
                          ridge = FALSE, gr.init=10, lambda.min=0.0001){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: gen_int_beta
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Generating the initial values using FMR or ridge fusion method.
  ## -----------------------------------------------------------------------------------------------------------------

  pp = p+q
  n = length(whole.data$data_y)
  sample.index  =  combn(n,2)
  a.matrix  =  kronecker(Matrix(t(apply(sample.index,2,dMatrixFun,n=n)),sparse = T),bdiag(diag(pp)))
  if(!ridge){
    fmr.bic = c()
    fmr.res0 = c()
    res.mle0.list = list()
    nC.vec = subgroup
    for (ii in 1:length(nC.vec)) {
      nC = nC.vec[ii]
      b00 = c(rep(-1,(nC*pp+nC)/2),rep(1,nC*pp+nC-floor((nC*pp+nC)/2)))
      res.mle  =  fmrs.mle(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
                           nComp =nC, disFamily = "norm",
                           # initCoeff = rnorm(nC*pp+nC),
                           initCoeff = b00,
                           initDispersion = rep(1,nC),
                           initmixProp = rep(1/nC, nC),nIterNR = 200)
      res.lam  =  fmrs.tunsel(y = whole.data$data_y, x =whole.data$data_x,delta = rep(0,length(whole.data$data_y)),
                              nComp =nC, disFamily = "norm",
                              initCoeff = c(coefficients(res.mle)),
                              initDispersion = dispersion(res.mle),
                              initmixProp = mixProp(res.mle),
                              penFamily = "adplasso",nIterNR = 200)
      res.var  =  fmrs.varsel(y = whole.data$data_y, x =whole.data$data_x, delta = rep(0,length(whole.data$data_y)),
                              nComp = ncomp(res.mle), disFamily = "norm",
                              initCoeff=c(coefficients(res.mle)),
                              initDispersion = dispersion(res.mle),
                              initmixProp = mixProp(res.mle),
                              penFamily = "adplasso",
                              lambPen = slot(res.lam, "lambPen"),nIterNR = 200)
      res.mle0 = res.var
      res.mle0.list[[ii]] = res.mle0

      fmr.res=apply(abs(res.mle0@residuals),1,which.min)
      fmr.resj=c()
      for (jj in 1:n) {
        fmr.resj[jj]=res.mle0@residuals[jj,fmr.res[jj]]
      }
      fmr.res0[ii] = sqrt(sum(fmr.resj^2))
      fmr.bic[ii] = res.mle0@BIC
    }
    fmr.coe = res.mle0@coefficients[-1,]
    fmr.w = apply(res.mle0@weights,1,which.max)
    fmr.r = apply(abs(res.mle0@residuals),1,which.min)

    # weight
    res.mle0  =  res.mle0.list[[which.max(fmr.bic)]]
    nC = nC.vec[which.max(fmr.bic)]
    beta.int.fmr = list()
    for (i in 1:n) {
      beta.int.fmr[[i]] = as.matrix(as.numeric(fmr.coe[,fmr.w[i]]))
    }
    beta.init=do.call(rbind,beta.int.fmr)
  }

  if(ridge){
    ridge.list  =  InitialBeta(whole.data, pp, n, lambda.min=lambda.min,
                               gr.init=gr.init, a.matrix=a.matrix)
    beta.init  =  ridge.list$beta.init
    beta.init.mat  =  t(matrix(beta.init,nrow = pp))
    nC  =  length(unique(apply(beta.init.mat, 1, function(a){floor(sum(a))})))
  }
  return(list(beta.init=beta.init,nC=nC))
}
