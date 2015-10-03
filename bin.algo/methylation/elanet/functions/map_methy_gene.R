map_methy_gene <- function(data,gene_set,probe2gene){
  
  all_gene = toupper(probe2gene[,2])
  gene_set = toupper(gene_set)
  
  ix = map(gene_set,all_gene)
  ix = ix[!is.na(ix)]
  
  if( length(ix)==0){
    return(NA)
  }else{
    probe_ix = match(probe2gene[ix,1],rownames(data))
    return(data[probe_ix,])
  }
  
}