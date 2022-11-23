#' Generate tuning parameters
#'
#' @author Mingyang Ren
#' @usage genelambda.obo(nlambda1=20,lambda1_max=0.5,lambda1_min=0.1,
#'                       nlambda2=5,lambda2_max=1.5,lambda2_min=0.1)
#' @description Generating a sequence of the tuning parameters (lambda1 and lambda2).
#'
#' @param nlambda1 The numbers of lambda 1.
#' @param lambda1_max The maximum values of lambda 1.
#' @param lambda1_min The minimum values of lambda 1.
#' @param nlambda2 The numbers of lambda 2.
#' @param lambda2_max The maximum values of lambda 2.
#' @param lambda2_min The minimum values of lambda 2.
#'
#' @return A sequence of the tuning parameters (lambda1, lambda2, and lambda3).
#' @export
#'
#' @examples
#' lambda <- genelambda.obo()
#' lambda

genelambda.obo = function(nlambda1=20,lambda1_max=0.5,lambda1_min=0.1,
                          nlambda2=5,lambda2_max=1.5,lambda2_min=0.1)
{

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: genelambda.obo
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Generating a sequence of the tuning parameters (lambda1 and lambda2).
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ nlambda1, nlambda2: The numbers of lambda 1 2.
  ## @ lambda1_min, lambda2_min: The minimum values of lambda 1 2.
  ## @ lambda1_max, lambda2_max: The maximum values of lambda 1 2.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ lambda: a sequence of the tuning parameters (lambda1 and lambda2).
  ## -----------------------------------------------------------------------------------------------------------------

  lambda1 = exp(seq(log(lambda1_max),log(lambda1_min),len= nlambda1))
  lambda2 =exp(seq(log(lambda2_max),log(lambda2_min),len= nlambda2))
  lambda = list(lambda1=lambda1,lambda2=lambda2)
  return(lambda)
}
