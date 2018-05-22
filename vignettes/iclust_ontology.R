C1 <- qusage::read.gmt("C:/RFactory/predsurv/vignettes/Gene sets/c1.all.v6.1.entrez.gmt")
gene_names <- readRDS("C:/RFactory/Rdata_brca/gene_names.RDS")

match("GATA3", gene_names$Hugo_Symbol)

gene_names$Entrez_Gene_Id[match("GATA3", gene_names$Hugo_Symbol)]

mg <- function(x){
  gene_names$Entrez_Gene_Id[match(x, gene_names$Hugo_Symbol)]
}


intClust1 <- unique(c(C1$chr17q23, mg("GATA3")  , mg("RPS6KB1") , mg("PPM1D,"), mg("PTRH2"), mg("APPBP2") ) )

intClust1 <-intClust1[!is.na(intClust1)]

intClust2 <- unique(c( C1$chr11q13, C1$chr11q14, mg("CCND1")  ,mg("EMSY")  ,mg("PAK1") ) )
intClust2 <-intClust2[!is.na(intClust2)]


intClustGO <- list(intClust1,intClust2 )
