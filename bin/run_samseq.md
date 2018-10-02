Run samseq for after/before seprately, and run for each cluster vs rest.

### Function for running SAMseq

``` r
run.SAMseq <- function(data,group1,group2,save.samfit=NULL){
   # data is matrix with expressions, samples as colnames, genes as rownames
   # OBS! has to be counts, samseq only takes integer values
   # groups, vectors with column indices in data for the groups to compare
   # save.samfit - optional to provide a file where the samfit object is saved!
  suppressMessages(library(samr))
  
   # change duplicate gene names
   rownames(data)<-change.duplicate.names(rownames(data))

   nSamp<-length(group1)+length(group2)
   y<-c(rep(1,length(group1)),rep(2,length(group2)))
   samfit<-SAMseq(data[,c(group1,group2)],y,resp.type = "Two class unpaired",fdr.output=0.05,genenames=rownames(data))
   if (!is.null(save.samfit)){ save(samfit,file=save.samfit) }

   nUp<-samfit$siggenes.table$ngenes.up
   nLo<-samfit$siggenes.table$ngenes.lo
   if (nUp==0 && nLo>0){
      ss<-samfit$siggenes.table$genes.lo
      type<-rep("down",nLo)
   }else if (nLo==0 && nUp>0){
      ss<-samfit$siggenes.table$genes.up
      type<-rep("down",nUp)
   }else if (nLo==0 && nUp==0){
      ss<-t(as.matrix(rep(NA,5)))
      colnames(ss)<-c("Gene ID", "Gene Name",  "Score(d)",  "Fold Change", "q-value(%)")
      type<-NA
   }else {
      ss <- rbind(samfit$siggenes.table$genes.up, samfit$siggenes.table$genes.lo)
      type<-c(rep("up",nUp),rep("down",nLo))
   }
   ss<-cbind(ss,type)
   return(ss)
}
```

Read data.
==========

``` r
# read rpkm values
C <-read.table("../data/ensembl_counts_filtered.csv", sep=",", header=T)

# read in meta data table
M <- read.table("../data/meta_data_with_clust.csv", sep=",",header=T)
```

Define clusters to run on.

``` r
after <- which(M$celltype == "graft")
before <- which(M$celltype == "Culture")

clusters.after <- split(after, as.character(M$MergeClust[after]))
clusters.before <- split(before, as.character(M$MergeClust[before]))
```

Run SAMseq for after
====================

``` r
for (cl in names(clusters.after)){
  outfile <- sprintf("../data/samseq/samseq_%s_vs_restAfter.txt",cl)
  if (!file.exists(outfile)){
    cells1 <- clusters.after[[cl]]
    cells2 <- setdiff(unlist(clusters.after), cells1)
    D<-run.SAMseq(C,cells1, cells2)
    write.table(D,file=outfile,sep="\t")
  }
}
```

Run SAMseq for before
=====================

``` r
for (cl in names(clusters.before)){
  outfile <- sprintf("../data/samseq/samseq_%s_vs_restBefore.txt",cl)
  if (!file.exists(outfile)){
    cells1 <- clusters.before[[cl]]
    cells2 <- setdiff(unlist(clusters.before), cells1)
    D<-run.SAMseq(C,cells1, cells2)
    write.table(D,file=outfile,sep="\t")
  }
}
```
