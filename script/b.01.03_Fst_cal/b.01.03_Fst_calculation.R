library("data.table")
geno <- fread("snp_calls.txt", header=FALSE)

names(geno) <- c("chr", "pos", "ref", "alt", "quality", "depth", paste0("l",1:17))


geno <- as.data.frame(geno)
for(i in 7:23){
  # replace slash and everything after it as nothing
  geno$newcol <- gsub("/.*", "", geno[,i] )
  # extract the line name
  nm <- names(geno)[i]
  # assign name for this allele
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a1")
  
  geno$newcol <- gsub(".*/", "", geno[,i] )
  names(geno)[ncol(geno)] <- paste0(nm, sep="_a2")
}

geno[geno == "."] <- NA
names(geno)

geno$p <- apply(geno[, 24:57], 1, function(x) {sum( x == 0) })
geno$p <- geno$p/34

geno$p1 <- apply(geno[, 24:41], 1, function(x) {sum(x ==0 )})
geno$p1 <- geno$p1/18

geno$p2 <- apply(geno[, 42:57], 1, function(x) {sum(x ==0) })
geno$p2 <- geno$p2/16

geno$fst <- with(geno, ((p1-p)^2 + (p2-p)^2)/(2*p*(1-p)) )
write.table(geno, "../fst.csv", sep=",", row.names = FALSE, quote=FALSE)
