# the rank 3 dsc

library(dscr)
library(flashr)
sessionInfo()

dsc_flashsim_r3=new_dsc("dsc-flashsim-r3","../output/dsc-flashsim-r3-files")

sim_r3_paper =function(){

  #set tau as in paper for defaults
  tau = 0.25
  n=150
  p=240

  l = matrix(0,nrow=n,ncol=3)
  f = matrix(0,nrow=p,ncol=3)
  l[1:10,1] = rnorm(10,0,2)
  l[11:60,2] = rnorm(50,0,1)
  l[61:150,3] = rnorm(90,0,0.5)

  f[1:80,1] = rnorm(80,0,.5)
  f[81:160,2] = rnorm(80,0,1)
  f[161:240,3] = rnorm(80,0,2)


  LF = l %*% t(f)

  Y = LF + sqrt(1/tau) * matrix(rnorm(n*p),nrow=n,ncol=p)
  return(list(Y=Y,LF=LF))
}

sim_r3_paper.wrapper = function(args){
  res = sim_r3_paper()
  return(list(meta = list(LF=res$LF),input = list(Y=res$Y)))
}


flash_r3.wrapper = function(input,args){
  f=flashr::flash_add_greedy(input$Y,Kmax=3)
  return(list(LFhat = flash_get_lf(f)))
}

flash_r3_pn.wrapper = function(input,args){
  f=flashr::flash_add_greedy(input$Y,Kmax=3,ebnm_fn=ebnm_pn)
  return(list(LFhat = flash_get_lf(f)))
}

udv = function(res){res$u %*% diag(res$d) %*% t(res$v)}
ssvd_r3.wrapper = function(input,args){
  res = ssvd::ssvd(input$Y,method = "theory",r=3)
  return(list(LFhat = udv(res)))
}

ssvd_r3.method.wrapper = function(input,args){
  res = ssvd::ssvd(input$Y,method = "method",r=3)
  return(list(LFhat = udv(res)))
}

svd_r3.wrapper = function(input,args){
  res = svd(input$Y,nu=3,nv=3)
  res$d = res$d[1:3]
  return(list(LFhat = udv(res)))
}

pmd_r3.wrapper = function(input,args){
  res = PMA::PMD(input$Y,center=FALSE,K=3)
  return(list(LFhat = udv(res)))
}

# pmd.cv.wrapper = function(input,args){
#   cv.out <- PMA::PMD.cv(input$Y, type="standard", sumabss=seq(0.1, 1, len=20),center=FALSE)
#   res <-  PMA::PMD(input$Y, type="standard", sumabs=cv.out$bestsumabs, K=1, v=cv.out$v.init,center=FALSE)
#   return(list(LFhat = res$d[1] * res$u[,1] %*% t(res$v[,1])))
# }


# this a helper function by Wei to do CV
CVPMD_softImpute=function(Y,c_s,K,fold = 10, method = "PMD"){
  N = dim(Y)[1]
  P = dim(Y)[2]
  colindex = matrix(sample(P,P),ncol = fold)
  rowindex = matrix(sample(N,N),ncol = fold)

  missing= array(0,dim = c(fold,N,P))
  foldindex = array(0,dim = c(fold,fold,2))
  for(i in 1:fold){
    for(j in 1:fold){
      foldindex[i,j,1] = i
      foldindex[i,j,2] = (i+j) %% fold
    }
  }
  foldindex[which(foldindex == 0)] = fold
  for(i in 1:fold){
    missing[i, , ] = Y
    for(j in 1:fold){
      missing[i,rowindex[,foldindex[j,i,1]],colindex[,foldindex[j,i,2]]] = NA
    }
    missing[i,,which(colSums(missing[i,,],na.rm = T) ==0)] = Y[,which(colSums(missing[i,,],na.rm = T) ==0)]
  }
  # c_s is the candicate of shrinkage parameter
  n_s = length(c_s)
  # rmse for each grids
  CVRMSE = rep(0,n_s)
  minrmse = Inf
  opt_s = 0
  # for each candidate, we run it N_sim times
  for(t_s in 1:n_s){
    # for each grid
    # each time we set the rmse to zeros
    rmse = rep(0,fold)
    for(i in 1:fold){
      if(method == "PMD"){
        res_log = capture.output({out = PMD(missing[i,,], sumabs = c_s[t_s], sumabsv = NULL, sumabsu = NULL,K = K)})
      }else{
        out = softImpute::softImpute(missing[i,,], rank.max = K,lambda = c_s[t_s])
      }
      if(length(out$d)==1){
        misshat = (out$d) * out$u %*% t(out$v)
      }else{
        misshat = out$u %*%  diag(out$d) %*% t(out$v)
      }
      for(j in 1:fold){
        # for each fold j
        rmse[i] = rmse[i] + sum((Y[rowindex[,foldindex[j,i,1]],colindex[,foldindex[j,i,2]]] -
                                   misshat[rowindex[,foldindex[j,i,1]],colindex[,foldindex[j,i,2]]])^2,na.rm = TRUE)
      }
    } #get the result for one run
    CVRMSE[t_s] = CVRMSE[t_s] + sqrt(sum(rmse)/(N*P))
    if(CVRMSE[t_s] < minrmse){
      minrmse = CVRMSE[t_s]
      opt_s = c_s[t_s]
    }
  }
  return(list(opt_s = opt_s, output = CVRMSE))
}



softImpute_r3.cv.wrapper = function(input,args){
  Y= input$Y
  K=3
  ngrids = args$ngrids
  fold = args$fold
  max = args$max

  N = dim(Y)[1]
  P = dim(Y)[2]
  c_s = seq(0,max,len=ngrids)
  cvout = CVPMD_softImpute(Y,c_s,K ,fold , method = "softImpute")
  res = softImpute::softImpute(Y, rank.max = K,lambda = cvout$opt_s)
  return(list(LFhat = res$d[1] * res$u %*% t(res$v)))
}



### Scores


rmse = function(a,b){return(sqrt(mean((a-b)^2)))}
rrmse.wrapper = function(data,output){
  return(list(rmse=rmse(data$meta$LF,output$LFhat), rrmse=rmse(data$meta$LF,output$LFhat)/sqrt(mean(data$meta$LF^2))))
}

add_scenario(dsc_flashsim_r3,"block-sparse-r3",sim_r3_paper.wrapper,seed=1:10)

add_method(dsc_flashsim_r3,"flash_r3",flash_r3.wrapper)
add_method(dsc_flashsim_r3,"flash_r3_pn",flash_r3_pn.wrapper)
add_method(dsc_flashsim_r3,"ssvd",ssvd_r3.wrapper)
add_method(dsc_flashsim_r3,"ssvd.method",ssvd_r3.method.wrapper)
add_method(dsc_flashsim_r3,"svd",svd_r3.wrapper)
add_method(dsc_flashsim_r3,"pmd",pmd_r3.wrapper)
# #add_method(dsc_flashsim,"pmd.cv",pmd.cv.wrapper)
add_method(dsc_flashsim_r3,"softimpute.cv",softImpute_r3.cv.wrapper,args=list(ngrids=10,fold=5,max=100))
add_method(dsc_flashsim_r3,"softimpute2.cv",softImpute_r3.cv.wrapper,args=list(ngrids=10,fold=5,max=10)) # finer grid near 0
#

add_score(dsc_flashsim_r3,rrmse.wrapper,"rmse")


res_r3=run_dsc(dsc_flashsim_r3)
save(res_r3,dsc_flashsim_r3,file="../output/dsc-flashsim-r3-files/res.RData")
