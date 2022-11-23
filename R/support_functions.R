dMatrixFun  =  function(indx, n){
  e.vec = matrix(0,n,1)
  e.vec[indx[1],] = 1
  e.vec[indx[2],] = (-1)
  return(e.vec)
}

mcp_d  =  function(x,lambda,a=3){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: mcp_d
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##            Calculating the derivative of the MCP
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input:
  ## @ x: a float value or a vector, the independent variable in the MCP.
  ## @ lambda: a float value, the tuning parameter in the MCP.
  ## @ a: a float value, regularization parameter in the MCP, the default setting is 3.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Output:
  ## @ rho: the derivative of the MCP.
  ## -----------------------------------------------------------------------------------------------------------------

  if(lambda!=0){
    rho  =  lambda*( 1 > abs(x)/( lambda*a ) )*( 1 - abs(x)/( lambda*a ))
  } else{
    rho=0
  }
  return(rho)
}

get.gr.num = function(size, beta.gfl,eta.gfl, q, sample.index, merge.all=F){

  ## -----------------------------------------------------------------------------------------------------------------
  ## The name of the function: get.gr.num
  ## -----------------------------------------------------------------------------------------------------------------
  ## Description:
  ##          compute the number of estimated subgroups.
  ## -----------------------------------------------------------------------------------------------------------------
  ## Input: It can refer to the output of the algorithm
  ## -----------------------------------------------------------------------------------------------------------------

  diff.gfl = apply(matrix(eta.gfl,nrow=q,ncol=size*(size-1)/2),2,sum)

  capgfl.matrix = matrix(0,nrow=size,ncol=size)
  if(length(which(diff.gfl==0))==0){
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    return(list(gr.num=size, capgfl.matrix2=capgfl.matrix2))
  }
  if(length(which(diff.gfl==0))==size*(size-1)/2){
    capgfl.matrix[upper.tri(capgfl.matrix)]=1
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    return(list(gr.num=1, capgfl.matrix2=capgfl.matrix2))
  }else{
    sample.index.gfl = sample.index[,which(diff.gfl==0)]
    if(length(which(diff.gfl==0))==1){
      capgfl.matrix[sample.index.gfl[1],sample.index.gfl[2]] = 1
    }else{
      for(i in 1:length(which(diff.gfl==0))){
        capgfl.matrix[sample.index.gfl[1,i],sample.index.gfl[2,i]] = 1
      }
    }
    capgfl.matrix2 = capgfl.matrix+t(capgfl.matrix); diag(capgfl.matrix2)=1
    group.num.gf = nrow(unique(capgfl.matrix2))

    if(merge.all){
      cap  =  capgfl.matrix2
      num_subgroup  =  unique(apply(cap, 1, function(a){which(a == 1)}))
      non_inter_list  =  list()
      vv  =  1
      non_inter  =  c(1:length(num_subgroup))
      repeat{
        a  =  num_subgroup[[non_inter[1]]]
        KK_k  =  setdiff(non_inter,non_inter[1])
        non_inter  =  c()
        i=1
        for (k2 in KK_k) {
          if(length(intersect(a,num_subgroup[[k2]])) > 0){
            a  =  union(a,num_subgroup[[k2]])
          } else {
            non_inter[i]  =  k2
            i=i+1
          }
        }
        non_inter_list[[vv]]  =  a
        vv  =  vv+1
        if(length(non_inter) == 0){break}
      }

      for (i in 1:dim(cap)[1]) {
        for (k in 1:length(non_inter_list)) {
          if(length(match(cap[i,],non_inter_list[[k]])) > 0){
            cap[i,non_inter_list[[k]]]  =  1
          }
        }
      }
      capgfl.matrix2  =  cap
      group.num.gf = nrow(unique(capgfl.matrix2))
    }

    return(list(gr.num=group.num.gf,capgfl.matrix2=capgfl.matrix2))
  }

}

