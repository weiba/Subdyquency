compartment_calculate=function(xl,compartment)
{
xl1=xl[intersect(row.names(xl),compartment[,1]),intersect(colnames(xl),compartment[,1])]
kk=weight_xl(xl1,compartment)
new_matrix=matrix(0,nrow=nrow(xl),ncol=ncol(xl),dimnames=list(row.names(xl),colnames(xl)))
new_matrix[intersect(row.names(new_matrix),row.names(kk)),intersect(colnames(new_matrix),colnames(kk))]=kk
return(new_matrix)
}
weight_xl=function(xl,compartment)
{
  relat_index=which(xl!=0,arr.ind = T)
  comp=compute_comp(compartment)
  for(i in 1:nrow(relat_index))
  {
    Gen1=row.names(xl)[relat_index[i,1]]
    Gen2=colnames(xl)[relat_index[i,2]]
    
      Gen1_com=compartment[which(Gen1==compartment[,1]),2]
      Gen2_com=compartment[which(Gen2==compartment[,1]),2]
    
   if(length(intersect(Gen1_com,Gen2_com))>0)
    {
      ll=intersect(Gen1_com,Gen2_com)
      xl[Gen1,Gen2]=max(as.numeric(comp[match(ll,comp[,1]),2]))
    }
    else
    {
      xl[Gen1,Gen2]=min(as.numeric(comp[,2]))
    }
  }
  return(xl)
  
}
compute_comp=function(compartment)
{
  comp=compartment[which(duplicated(compartment[,2])==FALSE),2]
  com=c()
  for(i in 1:length(comp))
  {
    com=rbind(com,cbind(as.character(comp[i]),length(which(compartment[,2]==comp[i]))))
  }
  com[,2]=as.numeric(com[,2])/max(as.numeric(com[,2]))
  
  return(com)
}
