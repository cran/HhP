#' Hierarchical Heterogeneity Regression Analysis.
#'
#' @author Mingyang Ren, Qingzhao Zhang, Sanguo Zhang, Tingyan Zhong, Jian Huang, Shuangge Ma. Maintainer: Mingyang Ren <renmingyang17@mails.ucas.ac.cn>.
#' @references Mingyang Ren, Qingzhao Zhang, Sanguo Zhang, Tingyan Zhong, Jian Huang, Shuangge Ma. 2022. Hierarchical Cancer Heterogeneity Analysis Based On Histopathological Imaging Features. Biometrics, <DOI: 10.1111/biom.13544>.
#' @usage HhP.reg(lambda, whole.data, n, q, p, beta.init,
#'                merge.all=FALSE, trace=FALSE, selection.sub=FALSE)
#'
#' @description The main function for Transfer learning for tensor graphical models.
#' @param lambda The sequences of the tuning parameters (lambda1 and lambda2).
#' @param whole.data The input data analyzed (a list including the response and design matrix).
#' @param n The sample size.
#' @param q The dimension of type 1 features.
#' @param p The dimension of type 2 features.
#' @param beta.init The Initial values of regression coefficients.
#' @param trace the logical variable, whether or not to output the number of identified subgroups during the search for parameters.
#' @param merge.all the logical variable, the default is F.
#' @param selection.sub the logical variable, the default is F.
#'
#'
#' @return A result list.
#' @export
#'
#' @import Matrix MASS fmrs
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

HhP.reg  =  function(lambda, whole.data, n, q, p, beta.init,
                     merge.all=FALSE, trace=FALSE, selection.sub=FALSE)
{

  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Searching and selecting the optional tuning parameters under
  ##            the adaptive BIC-type criterion using the proposed method.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ lambda: The sequences of the tuning parameters (lambda1 and lambda2).
  ## @ whole.data: The input data analyzed (a list including the response and design matrix).
  ## @ n: The sample size.
  ## @ q: The dimension of type 1 features.
  ## @ p: The dimension of type 2 features.
  ## @ beta.init: The Initial values of regression coefficients.
  ## @ trace: the logical variable, whether or not to output the number of identified subgroups during the search for parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------
  ## Output:
  ## A list including results corresponding all choices of given tuning parameters.
  ## ------------------------------------------------------------------------------------------------------------------------------------------


  pp  =  q+p
  sample.index  =  combn(n,2)
  a.matrix  =  kronecker(Matrix(t(apply(sample.index,2,dMatrixFun,n=n)),sparse = T),bdiag(diag(pp)))
  one.matrix  =  bdiag(rep(list(rep(1,pp)),(n*(n-1)/2)))
  one.matrix.main  =  bdiag(rep(list(rep(1,q)),(n*(n-1)/2)))
  one.matrix.GE  =  bdiag(rep(list(rep(1,p)),(n*(n-1)/2)))


  lambda1 = lambda$lambda1
  lambda2 = lambda$lambda2
  L1 = length(lambda1)
  L2 = length(lambda2)
  L = L1+L2

  lam2 = lambda2
  lam1 = lam1=rep(0,length(lam2))
  admmres2=vector(mode="list",length=length(lam2))
  bic.var2 = c()
  gr.nn2  =  c()
  for(l2 in 1:L2){
    if(trace){if(l2%%5==1){cat('-----------',l2,'-th lambda--------------\n')}}
    ad = ADMM(n, q, p, whole.data, beta.init, lam1[l2], lam2[l2], merge.all=merge.all,
              sample.index=sample.index,a.matrix=a.matrix, one.matrix=one.matrix,
              one.matrix.main=one.matrix.main, one.matrix.GE=one.matrix.GE)
    bic.var2[l2] = ad$BIC.var
    RI.main=ad$num.main
    RI.GE=ad$num.GE
    admmres2[[l2]]  =  ad
    gr.nn2[l2]  =  RI.main$gr.num
    if(trace){
      print(c(RI.main$gr.num,RI.GE$gr.num))
    }
  }

  if(sum(gr.nn2 > 1) > 0){
    abic.n2 = which(bic.var2 == min(bic.var2[(!is.na(bic.var2)) & (gr.nn2 > 1)]))[1]
  } else {
    abic.n2 = which(bic.var2 == min(bic.var2[(!is.na(bic.var2)) & (gr.nn2 > 0)]))[1]
  }

  lam20=lam2[abic.n2][1]
  lam1=lambda1
  lam2=rep(lam20,length(lam1))
  admmres=vector(mode="list",length=length(lam1))
  bic.var = c()
  for(l1 in 1:L1){
    if(trace){if((L2+l1)%%5==1){cat('-----------',L2+l1,'-th lambda--------------\n')}}
    ad = ADMM(n, q, p, whole.data, beta.init, lam1[l1], lam2[l1], merge.all=merge.all,
              sample.index=sample.index,a.matrix=a.matrix, one.matrix=one.matrix,
              one.matrix.main=one.matrix.main, one.matrix.GE=one.matrix.GE)
    bic.var[l1] = ad$BIC.var
    RI.main=ad$num.main
    RI.GE=ad$num.GE
    admmres[[l1]]  =  ad
    if(trace){
      print(c(RI.main$gr.num,RI.GE$gr.num))
    }
  }
  abic.n = which(bic.var == min(bic.var[!is.na(bic.var)]))[1]
  beta.over = admmres[[abic.n]]$beta.gfl
  RI.main  =  admmres[[abic.n]]$RI.main
  RI.GE  =  admmres[[abic.n]]$RI.GE
  result = list(beta.over=beta.over,bic.var=bic.var,abic.n=abic.n, admmres=admmres,
                bic.var2=bic.var2, abic.n2=abic.n2, admmres2=admmres2)

  if(selection.sub){
    admm.sum  =  admmres
    num.sum  =  matrix(0,length(admm.sum),2)
    BIC  =  c()
    for (j in 1:length(admm.sum)) {
      admmj  =  admm.sum[[j]]
      BIC[j]  =  admmj$BIC.var
      num.sum[j,1]  =  admmj$num.main$gr.num
      num.sum[j,2]  =  admmj$num.GE$gr.num
    }
    select.num  =  which(num.sum[,1] > 1 & num.sum[,2] > num.sum[,1])
    if(length(select.num) == 0){
      result = list(beta.over=beta.over,bic.var=bic.var,abic.n=abic.n, admmres=admmres,
                    bic.var2=bic.var2, abic.n2=abic.n2, admmres2=admmres2)
    }
    if(length(select.num) > 0){
      abic.n  =  which(BIC == min(BIC[select.num]))
      beta.over = admmres[[abic.n]]$beta.gfl
      RI.main  =  admmres[[abic.n]]$RI.main
      RI.GE  =  admmres[[abic.n]]$RI.GE
      result = list(beta.over=beta.over,bic.var=bic.var,abic.n=abic.n, admmres=admmres,
                    bic.var2=bic.var2, abic.n2=abic.n2, admmres2=admmres2)
    }
  }
  return(result)
}


