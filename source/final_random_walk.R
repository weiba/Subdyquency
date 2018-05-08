library(igraph)
main_function=function(patMutMatrix,patOutMatrix,influenceGraph,a=0.5)
{pre_result=preprocess_func(patMutMatrix,patOutMatrix,influenceGraph)
xl=.buildBipartiteGraph_1(pre_result$patMut,pre_result$patOut,pre_result$infGraph)
di=compute_ratio(pre_result$patMut,pre_result$patOut)
new_matrix=compartment_calculate(xl,compartment)
patMutMatrix=patMutMatrix[,colnames(new_matrix)]
mut_rank=patient_mutandout(patMutMatrix,di$mut_ratio,di$out_ratio,new_matrix,a)
return(mut_rank)
}
preprocess_func=function(patMutMatrix,patOutMatrix,influenceGraph)
{
  diag(influenceGraph)=0
  if(length(which(colSums(patMutMatrix)==0))>0)
  {
    patMutMatrix=patMutMatrix[,-which(colSums(patMutMatrix)==0)]
  }
  if(length(which(colSums(patOutMatrix)==0))>0)
  {
    patOutMatrix=patOutMatrix[,-which(colSums(patOutMatrix)==0)]
  }
  OutInter=intersect(row.names(influenceGraph),colnames(patOutMatrix))
  MutInter=intersect(colnames(influenceGraph),colnames(patMutMatrix))
  patMutMatrix=patMutMatrix[,MutInter]
  patOutMatrix=patOutMatrix[,OutInter]
  pats=intersect(row.names(patMutMatrix),row.names(patOutMatrix))
  patMutMatrix=patMutMatrix[pats,]
  patOutMatrix=patOutMatrix[pats,]
  influenceGraph=influenceGraph[OutInter,MutInter]
  list(patMut=patMutMatrix,patOut=patOutMatrix,infGraph=influenceGraph)
}
.buildBipartiteGraph_1=function(patMut,patOut,infG)
{
  BiG=matrix(0,nrow = ncol(patOut),ncol = ncol(patMut))
  row.names(BiG)=colnames(patOut)
  colnames(BiG)=colnames(patMut)
  
  Candidate_Driver=colnames(patMut)
  for(i in 1:length(Candidate_Driver))
  {
    RelateOut=row.names(infG)[which(infG[,Candidate_Driver[i]]==1)]
    RelatePat=which(patMut[,Candidate_Driver[i]]==1)
    if(length(RelatePat)==0)
    {next}
    for(j in 1:length(RelatePat))
    {
      pat_events=intersect(RelateOut,colnames(patOut)[which(patOut[RelatePat[j],]==1)])
      if(length(pat_events)>0)
      {BiG[pat_events,Candidate_Driver[i]]=1
      #message('daozhelila')}
      }
    }
  }
  count_mutated <- rowSums(BiG)
  count_outlier <- colSums(BiG)
  del_row <- which(count_mutated == 0)
  del_col <- which(count_outlier == 0)
  if (length(del_row) > 0) {
    BiG <- BiG[-del_row, ]
  }
  if (length(del_col) > 0) {
    BiG <- BiG[, -del_col]
  }    
  BiG
}

compute_ratio=function(patMutMatrix,patOutMatrix)
{mut_ratio=c()
for(i in 1:ncol(patMutMatrix))
{
  mut_ratio_value=length(which(patMutMatrix[,i]==1))/nrow(patMutMatrix)
  mut_ratio=rbind(mut_ratio,cbind(Mut_Gen=colnames(patMutMatrix)[i],mut_ratio_value))
}

out_ratio=c()
out_ratio=cbind(colnames(patOutMatrix),1/ncol(patOutMatrix))
#out_ratio=cbind(colnames(patOutMatrix),1/1000)
intt=intersect(out_ratio[,1],mut_ratio[,1])
out_va=mut_ratio[match(intt,mut_ratio[,1]),]
out_ratio[match(intt,out_ratio[,1]),2]=out_va[,2]

list(mut_ratio=mut_ratio,out_ratio=out_ratio)
}


patient_mutandout=function(patMutMatrix,mut_ratio,out_ratio,influenceGraph,a)
{
  mut_matrix=t(patMutMatrix)#####row is gene and col is patients
  #mut_mut=influenceGraph[colnames(patMutMatrix),colnames(patMutMatrix)]
  out_matrix=matrix(0,nrow = nrow(influenceGraph),ncol = nrow(patMutMatrix),dimnames = list(row.names(influenceGraph),row.names(patMutMatrix)))
  for(i in 1:nrow(mut_matrix))
  {
    mutname=row.names(mut_matrix)[i]
    pats=which(mut_matrix[i,]!=0)
    m_ratio=as.numeric(mut_ratio[which(mutname==mut_ratio[,1]),2])
   # browser()
    if(length(pats)>0)
    {mut_matrix[i,pats]=m_ratio
    }
    else
    {
      next
    }
  }
  
  for(p in 1:nrow(patMutMatrix))
  {
    mutnames=colnames(patMutMatrix)[which(patMutMatrix[p,]!=0)]
    #     if(length(mutnames)==0)
    #     {
    #       next
    #     }
    for(i in 1:length(mutnames))
    {
      if(length(which(mutnames[i]==colnames(influenceGraph)))>0)
      {outname=row.names(influenceGraph)[which(influenceGraph[,mutnames[i]]==1)]
      out_matrix[outname,p]=as.numeric(out_ratio[match(outname,out_ratio),2])
      
      }
    }
  }
  
  rowsum=rowSums(influenceGraph)
  colsum=colSums(influenceGraph)
  influenceGraph_row=influenceGraph
  influenceGraph_col=influenceGraph
  for(i in 1:length(rowsum))
  {
    if(rowsum[i]!=0)
    {influenceGraph_row[i,]=influenceGraph[i,]/rowsum[i]}
  }
  for(j in 1:length(colsum))
  {
    if(colsum[j]!=0)
    {influenceGraph_col[,j]=influenceGraph[,j]/colsum[j]}
  }

  m_1=a*mut_matrix+(1-a)*t(influenceGraph_row)%*%out_matrix
   o_1=a*out_matrix+(1-a)*influenceGraph_col%*%m_1
   m_f=a*mut_matrix+(1-a)*t(influenceGraph_row)%*%o_1
 finalScore=rowSums(m_f)
  message(a)
  #finalScore=rowSums(m_1)
  finalScore=finalScore[order(finalScore,decreasing = T)]
  return(finalScore)
}
