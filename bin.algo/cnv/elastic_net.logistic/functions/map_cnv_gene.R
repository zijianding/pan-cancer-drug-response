map_cnv_gene <- function(data,gene_set){
  
  all_gene = rownames(data)
  genes = sapply(all_gene,function(x){tmp=strsplit(x,split="\\|"); tmp=unlist(tmp);return(tmp[1])})
  
  genes = toupper(genes)
  gene_set = toupper(gene_set)
  
  ix = match(gene_set,genes)  
  ix = ix[!is.na(ix)]
  
  if( length(ix)==0){
    return(NA)
  }else{
    return(data[ix,])
  }
  
}