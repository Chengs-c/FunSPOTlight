ct_topic_fun=function(U,cell.type.factor,ct_mode){
  
  U_topic=t(U)
  cell.type.factor.all=unique(cell.type.factor)
  ct_topic_fun=matrix(nrow = nrow(U_topic),ncol = length(cell.type.factor.all),0)
  colnames(ct_topic_fun)=cell.type.factor.all
  if(ct_mode=="median"){
  for(i in 1:(nrow(U_topic))){
    for (j in 1:length(cell.type.factor.all)) {
      ct_topic_fun[i,j]=median(U_topic[i,cell.type.factor==cell.type.factor.all[j]])
    }
  }
  }else{
    for(i in 1:(nrow(U_topic))){
      for (j in 1:length(cell.type.factor.all)) {
        ct_topic_fun[i,j]=mean(U_topic[i,cell.type.factor==cell.type.factor.all[j]])
      }
    }
  }
  colnames(ct_topic_fun) <- gsub("[[:punct:]]|[[:blank:]]",
                                        ".",
                                        colnames(ct_topic_fun))
  
  return(ct_topic_fun)
}
