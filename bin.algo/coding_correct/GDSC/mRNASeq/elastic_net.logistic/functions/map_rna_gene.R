#data as a data frame
#each column as a sample, as col.name
#each row as a feature, as row.name
#RNA-seq row.name: "TP53|ID", etc.
#gene_set: a character vector


map_rna_gene <- function(data,gene_set){
  
  all_gene = rownames(data)
  all_gene =  strsplit( x=all_gene, split="\\|"  ) 
  all_gene = lapply(all_gene,function(x) return(x[[1]][1]) )
  all_gene.symbol = toupper(all_gene.symbol)
  
  gene_set = toupper(gene_set)
  ix = match(gene_set,all_gene.symbol)
  ix = ix[!is.na(ix)]
  
  if(length(ix)==0)
  {
    print("No genes were mapped!")
    return(NA)
  }
  if(length(ix)>0)
  {
    cat(length(ix),"genes are mapped!",sep=" ")
    
    data = data[ix,]
    return(data)
    
  }
  
  
}