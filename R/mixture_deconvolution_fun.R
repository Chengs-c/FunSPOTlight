mixture_deconvolution_fun=function(ST.U,ct_topic_fun,cell.type.factor,min_count){
  
  ST.U=t(ST.U)
  cell.type.factor.all=unique(cell.type.factor)
  decon_mtrx=matrix(nrow = length(cell.type.factor.all), 
                      ncol = ncol(ST.U),
                      0)
  rownames(decon_mtrx)=colnames(ct_topic_fun)
  
  #print("Deconvoluting spot to")
  #total=ncol(profile_mtrx)
  #pb=txtProgressBar(min=0,max=total,style = 3)
  
  for (i in 1:ncol(ST.U)) {
    ##nnls得到细胞类型组成
    nnls_pred=nnls::nnls(A=ct_topic_fun,b=ST.U[,i])
    weights=nnls_pred$x
    
    #得到比例
    comp=weights/sum(weights)
    
    #去除最小
    comp[comp<min_count]=0
    weights[comp<min_count]=0
    
    comp_prop=comp/sum(comp)
    comp_prop[is.na(comp_prop)]=0
    
    decon_mtrx[,i]=comp_prop
    #setTxtProgressBar(pb,i)
  
  }
  
  #close(pb)
  
  return(decon_mtrx)
}
