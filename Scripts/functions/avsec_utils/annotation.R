##' Convert hgnc symbol (Human - org.Hs.eg.db) into entrez id. In case a single symbol has multiple entrez_id's use the first one.
##' 
##' @param genes Character vector of gene symbols to be converted
##' @return Character vector of entrez id's. It's NA in case it wasn't found.
symbol2entrez_id <- function(genes){
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  dt <- select(org.Hs.eg.db, keys=genes, columns=c("SYMBOL","ENTREZID"), keytype="SYMBOL") %>% as.data.table
  dt <- dt[, .(ENTREZID = ENTREZID[1]), by = SYMBOL]
  return_genes <- dt$ENTREZID

  if(!any(duplicated(genes)) && any(duplicated(return_genes))) warning("input genes are unique. Output genes are not!")
  
  return(return_genes)
}


##' Get the annotation hash table for a hgnc (gene name) symbol in Human (org.Hs.eg.db)
##'
##' @param genes Character vector of hgnc gene names.
##' @param features Feature vector of possible features. Check with `columns(org.Hs.eg.db)` all the available features.
##' @param keytype In which format are the given genes?
##' 
##' Details:
##' - Rows containing any NA's are omitted
##' - used database: org.Hs.eg.db
gene_symbol_anno_table <- function(genes, features = c("SYMBOL", "ENSEMBL","ENTREZID"), keytype="SYMBOL") {
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  # select is overriden by other packages!!!
  dt <- AnnotationDbi::select(org.Hs.eg.db, keys=unique(genes), columns=features, keytype=keytype) %>% as.data.table
  dt <- na.omit(dt)

  ## number of omitted genes
  gene_omitted <- sum(!genes %in% dt[, keytype, with = F][[1]])
  if(gene_omitted > 0) message(gene_omitted, " genes with no annotation were omitted")
  
  setkeyv(dt, keytype)
  return(dt[, .SD[1], by = keytype])
}


## convert ensembl id 2 symbol
ensembl2symbol <- function(genes) {
  dt <- gene_symbol_anno_table(genes, features = c("SYMBOL", "ENSEMBL"), keytype = "ENSEMBL")
  ## stopifnot(nrow(dt) == length(genes), all(dt[, SYMBOL ==genes]))
  setkey(dt, "ENSEMBL")
  dt <- dt[genes]
  return_genes <- dt[, SYMBOL]
  if(!any(duplicated(genes)) && any(duplicated(return_genes))) warning("input genes are unique. Output genes are not!")
  return(return_genes)
}

## convert ensembl id 2 symbol
symbol2ensembl <- function(genes) {
  dt <- gene_symbol_anno_table(genes, features = c("SYMBOL", "ENSEMBL"), keytype = "SYMBOL")
  ## stopifnot(nrow(dt) == length(genes), all(dt[, SYMBOL ==genes]))
  setkey(dt, "SYMBOL")
  dt <- dt[genes]
  return(dt[, ENSEMBL])
}

CODONS <- c("TTT", "TTC", "TTA", "TTG", "TCT", "TCC", "TCA", "TCG", "TAT", 
            "TAC", "TGT", "TGC", "TGG", "CTT", "CTC", "CTA", "CTG", "CCT", 
            "CCC", "CCA", "CCG", "CAT", "CAC", "CAA", "CAG", "CGT", "CGC", 
            "CGA", "CGG", "ATT", "ATC", "ATA", "ATG", "ACT", "ACC", "ACA", 
            "ACG", "AAT", "AAC", "AAA", "AAG", "AGT", "AGC", "AGA", "AGG", 
            "GTT", "GTC", "GTA", "GTG", "GCT", "GCC", "GCA", "GCG", "GAT", 
            "GAC", "GAA", "GAG", "GGT", "GGC", "GGA", "GGG")    
