fpca_single_cell1=function(geno,t,nbasis,genobasis,cell.type.factor,nfpca){
  cell.type.factor.all=unique(cell.type.factor)
  Phi <- eval.basis(t, genobasis)
  geno.max = max(geno)
  # geno = geno / geno.max
  #Phi的维度为特征数量（基因数量）*基函数数量
  #geno的维度为样本数*特征数量（基因数量）
  geno <- t(geno)
  U <- geno %*% Phi %*% ginv(t(Phi) %*% Phi)
  fd_obj1=fda::fd(coef = t(U),basisobj = genobasis)
  #C=U
  fd_obj1_eval=eval.fd(t,fd_obj1)
  fpca1=pca.fd(fd_obj1,nharm = nfpca)
  #fpca后的基函数
  fpca1_basis_eval=Phi%*%fpca1$harmonics$coefs
  
  fpca_matrix=fpca1_basis_eval
  # for (i in 1:length(cell.type.factor.all)) {
  #   U_celltype=U[which(cell.type.factor==cell.type.factor.all[i]),]
  #   fd_celltype_obj=fd(coef = t(U_celltype),basisobj = genobasis)
  #   fd_celltype_eval=eval.fd(t,fd_celltype_obj)
  #   fpca_celltype=pca.fd(fd_celltype_obj,nharm=nbasis)
  #   fpca_celltype_eval=Phi%*%fpca_celltype$harmonics$coefs[,(nbasis-1):nbasis]
  #   fpca_matrix=cbind(fpca_matrix,fpca_celltype_eval)
  # }
  # 
  # for (i in 1:length(cell.type.factor.all)) {
  #   U_celltype=U[which(cell.type.factor==cell.type.factor.all[i]),]
  #   fd_celltype_obj=fd(coef = t(U_celltype),basisobj = genobasis)
  #   fd_celltype_eval=eval.fd(t,fd_celltype_obj)
  #   fpca_celltype=pca.fd(fd_celltype_obj,nharm=nbasis)
  #   fpca_celltype_eval=Phi%*%fpca_celltype$harmonics$coefs[,1]
  #   fpca_matrix=cbind(fpca_matrix,fpca_celltype_eval)
  # }
  # 
  return(fpca_matrix)
  
}
