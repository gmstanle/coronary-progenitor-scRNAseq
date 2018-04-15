library(data.table)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(plyr)
library(foreach)
library(doMC)
library(ggplot2)
library(Matrix)
library(dplyr)
library(ggrepel)
library(cowplot)
library(Rmisc)
library(Seurat)
library(parallel)

## Modified Seurat functions
SubsampByAveraging <- function(data, ident.use, max.cells.per.ident = 10){
  for(id in ident.use) {
    cells.id <- data@cell.names[data@ident==id]
    if(length(cells.id) > max.cells.per.ident){
      id.averages <- data.frame(cell.name = paste0(id, 1:max.cells.per.ident))
      n.cells.avg <- ceiling(length(cells.id) / max.cells.per.ident)

      cells.used <- character(0)
      for(i in 1:n.cells.out){
        pool <- cells.id[!cells.id %in% cells.used]
        if(length(pool) >= n.cells.avg) {
          cells.avg <- sample(pool, n.cells.avg)

        }
      }
    }
  }
}

# Get the mean correlation of all genes to a small subset (for example, a set of known IEGs)
# data is a Seurat object
GetGenesByCorrelation <- function(seur, target.genes, threshold.cor = 0.2, genes.use = NA){
  gene.score <- GetGeneGeneCorrelations(seur, target.genes, genes.use)
  print(plot(gene.score[!names(gene.score) %in% target.genes][1:100]))
  return(unique(c(target.genes, names(gene.score[gene.score > threshold.cor]))))
}

GetGeneGeneCorrelations <- function(seur, target.genes, genes.use=NA){
  mat.log10 <- as.matrix(t(seur@scale.data))
  if(is.na(genes.use)) genes.use <- colnames(mat.log10)

  genes.use <- genes.use[genes.use %in% colnames(mat.log10)[which(apply(mat.log10, 2, sd)>0)]]

  gene.cors <- sapply(target.genes, function(x) as.numeric(cor(mat.log10[, x], mat.log10[, genes.use])))
  rownames(gene.cors) <- genes.use
  gene.cors <- gene.cors[!is.na(gene.cors[, 1]), ]
  gene.score <- rowMeans(gene.cors)
  gene.score <- sort(gene.score, decreasing = T)
  return(gene.score)
}

# Calculate rPCA with PC loadings recalculated as sum of top 30 + minus top 30 - genes by loading
CalcRPCA <- function(data, genes.use = NA, ncp=10, num.sig.genes = 30, subsample.genes = NA){
  if(is.na(genes.use)) genes.use <- rownames(data@data)
  data.pca <- t(as.matrix(data@data))

  if(!is.na(subsample.genes)) data.pca <- subsample.numgenes(data.pca, ng.min = 2200)

  data.pca <- data.pca[, genes.use]
  # using non-normalized data seems to work better for rPCA
  pcalist <- do.robpca(data.pca,ncp = ncp, from.mat = T)
  loadings <- pcalist[[3]]

  # Scale each gene to max of 1
  data.scaled <- data.pca / colMaxs(data.pca)

  # using Seurat-scaled gene expr values does not seem to work as well, probably due to outlier cells/groups
  # data.scaled <- data.seuratNorm / colMaxs(data.seuratNorm )
  library(pbapply)
  dims <- paste0('PC', 1:ncp)
  pc.sigs <- as.data.frame(do.call(cbind, lapply(dims, function(dim) {
    setorderv(loadings,dim,-1)
    genes.pos <- loadings[,gene][1:num.sig.genes]
    setorderv(loadings,dim,1)
    genes.neg <- loadings[,gene][num.sig.genes:1]

    # generate signature value of top PC +/- correlated genes
    # normalize every gene
    # Renames colums to PCx.pos, PCx.neg, and PCx for the xth PC
    cast.sig <- as.data.frame(cbind(rowSums(data.scaled[, genes.pos]), rowSums(data.scaled[, genes.neg])))
    cast.sig$pc.score <- cast.sig[, 1] - cast.sig[, 2]
    colnames(cast.sig) <- c(paste0(dim,'.pos'), paste0(dim,'.neg'), dim)

    return(as.matrix(cast.sig))
  }
  )))

  # jam the new PC signatures into the Seurat object
  pc.sigs <- pc.sigs[match(data@cell.names, rownames(pc.sigs)), ]
  pc.sigs <- as.matrix(pc.sigs[, grepl("PC[0-9]+$", colnames(pc.sigs))]) #only works for < 100 PCs
  data <- SetDimReduction(data, 'pca', 'cell.embeddings', new.data = pc.sigs)
}

## End of modified Seurat functions

# useful for Seurat objects (use: group_averages(tiss@data, tiss@ident))
group_averages <- function(mat, groups){
  group_names = unique(groups)
  means = matrix(0, dim(mat)[1], length(group_names))
  colnames(means) = group_names
  rownames(means) = rownames(mat)
  for(group in group_names){
    means[,group] = Matrix::rowMeans(mat[,groups == group,drop=FALSE])
  }
  means
}


norm.log.counts <- function(counts) {
  norm.fact <- colSums(counts)
  counts.norm <- t( apply(counts, 1, function(x) {x/norm.fact*1000000+1}))
  counts.log <- log2(counts.norm)
  counts.log
}

load.count.data <- function(fname, mincount = 1e5, min.cells=0, min.num.genes = 200, omit.bad.cells=TRUE) {
  counts <- fread(fname)

  cell.names <- counts[['V1']]
  length(cell.names)

  counts <- as.matrix(counts[,-1,with=F])
  rownames(counts) <- cell.names
  counts <- t(counts)

  gene.names = rownames(counts)
  ercc.counts <- counts[grepl("^ERCC-", gene.names),]
  tech.counts <- counts[grepl("^N_", gene.names),]
  ercc.names = rownames(ercc.counts)
  tech.names = rownames(tech.counts)
  counts <- counts[!rownames(counts) %in% c(ercc.names, tech.names),]
  counts[is.na(counts)] <- 0


  # omit genes detected in less than min.cells cells
  if(min.cells > 0) counts <- counts[rowSums(counts) >= min.cells,]

  # remove low-quality cells
  plot(log10(colSums(counts) + 1), colSums(counts > 0), breaks=100)
  abline(v=mincount, col="red")
  abline(h=min.num.genes, col = 'red')

  if(omit.bad.cells){
    ercc.counts <- ercc.counts[,colSums(counts)>mincount];
    last3.counts <- last3.counts[,colSums(counts)>mincount];
    counts <- counts[,colSums(counts)>mincount];
  }

  list(counts=counts, ercc.counts=ercc.counts, tech.counts=tech.counts)
}

# qc.stats <- function(){}

# returns dataframe with gene expression and cell info, from matrix of expression and
# data.table of cell info
# return: cell info is in first N columns, columns thereafter are gene expression
# returns NA for gene expr of cells that are in info.dt but not in expr.mat
get.expr.cellinfo <- function(expr.mat, info.dt, info.cols, cells.use, to.df = F){

  expr.mat = expr.mat[cells.use, , drop = F]
  info.dt = info.dt[cell.name %in% cells.use]
  dt.return = cbind(info.dt[, c('cell.name', info.cols), with = F],
                    expr.mat[match(info.dt[['cell.name']], rownames(expr.mat)), ])

  colnames(dt.return)[(length(info.cols) + 2):ncol(dt.return)] <- colnames(expr.mat)

  if(to.df){
    dt.return = as.data.frame(dt.return)
    rownames(dt.return) <- dt.return$cell.name
    dt.return <- dt.return[, 2:ncol(dt.return)]
  }
  return(dt.return)
}

seqplate.plt <- function(g.plot, counts, libinfo, cells.use, log.trans = T, title = '', range = NA){
  library(rasterVis)
  library(lattice)
  rowcontam.tocast = get.expr.cellinfo(counts[, g.plot, drop = F],
                                       libinfo,
                                       c('i5.index','i7.index', 'seq.plate'),
                                       cells.use = cells.use)
  plate.expr=as.data.frame(dcast.data.table(rowcontam.tocast, i5.index ~ i7.index, value.var = g.plot, fill = NA))
  rownames(plate.expr)=plate.expr[,'i5.index']
  plate.expr=as.matrix(plate.expr[,2:ncol(plate.expr)])
  if(log.trans) plate.expr = log10(plate.expr + 1)
  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  myTheme <- BTCTheme()
  myTheme$panel.background$col = 'grey30'
  if(is.na(range)){
    p = levelplot(t(plate.expr),col.regions=color.palette(100),scales=list(x=list(rot=45)),
                  main=title, par.settings = myTheme)
  } else {
    p = levelplot(t(plate.expr),col.regions=color.palette(100),scales=list(x=list(rot=45)),
                  main=title, par.settings = myTheme, at=seq(range[1], range[2], length.out=100))
  }
  return(p)
}

# plot some metadata variable mapped to a sequencing plate (seq.row by seq.col)
# variable must be found in the libinfo data.table
# empty wells will be NA; NA values will be colored grey
# If there are redundant samples (e.g. if >1 seq plate is included) then you will get warning
# "Aggregate function missing, defaulting to 'length"
seqplate.plt.meta <- function(var.plot, libinfo, log.trans = T, title = '', range = NA){
  library(rasterVis)
  library(lattice)

  plate.expr=as.data.frame(dcast.data.table(libinfo, seq.row ~ seq.col, value.var = var.plot, fill = NA))
  rownames(plate.expr)=plate.expr[, 'seq.row']
  plate.expr=as.matrix(plate.expr[, 2:ncol(plate.expr)])
  if(log.trans) plate.expr = log10(plate.expr + 1)

  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  myTheme <- BTCTheme()
  myTheme$panel.background$col = 'grey40'
  if(length(range) < 2){
    p = levelplot(t(plate.expr),col.regions=color.palette(100),scales=list(x=list(rot=45)),
                  main=title, par.settings = myTheme)
  } else {
    p = levelplot(t(plate.expr),col.regions=color.palette(100),scales=list(x=list(rot=45)),
                  main=title, par.settings = myTheme, at=seq(range[1], range[2], length.out=100))
  }

  p
}

calc.mode.fail.dist <- function(data.counts, n.cores = 8, libloc, ...) {
  load.scde(libloc)
  # data.counts is a melted data.table with a column for cell names and a column for genes and a column for counts
  gene.counts <- cast.for.scde(data.counts)

  # gene.counts is a data.frame with genes in teh rows and cells in the columns. you can just start here and rewrite the function to bypass data.counts.
  o.ifm <- scde.error.models(counts=gene.counts,n.cores=n.cores,threshold.segmentation=T,save.crossfit.plots=F,save.model.plots=F,verbose=1, ...)
  valid.cells <- o.ifm$corr.a >0
  o.prior <- scde.expression.prior(models=o.ifm,counts=gene.counts,length.out=400,show.plot=T)
  o.fpm <- scde.expression.magnitude(o.ifm,counts=gene.counts);

  o.fail.curves <- scde.failure.probability(o.ifm,magnitudes=log((10^o.prior$x)-1))
  par(mfrow=c(1,1),mar = c(3.5,3.5,0.5,0.5), mgp = c(2.0,0.65,0), cex = 1);
  plot(c(),c(),xlim=range(o.prior$x),ylim=c(0,1),xlab="expression magnitude (log10)",ylab="drop-out probability")
  invisible(apply(o.fail.curves,2,function(y) lines(x=o.prior$x,y=y,col="orange")))

  p.self.fail <- scde.failure.probability(models=o.ifm,counts=gene.counts)
  cell.names <- colnames(gene.counts); names(cell.names) <- cell.names;

  # reclculate posteriors with the individual posterior modes
  jp <- scde.posteriors(models=o.ifm,gene.counts,o.prior,return.individual.posterior.modes=T,n.cores=n.cores)
  # find joint posterior modes for each gene - a measure of MLE of group-average expression
  jp$jp.modes <- log(as.numeric(colnames(jp$jp)))[max.col(jp$jp)]
  p.mode.fail <- scde.failure.probability(models=o.ifm,magnitudes=jp$jp.modes)
  # weight matrix
  matw <- 1-sqrt(p.self.fail*sqrt(p.self.fail*p.mode.fail))
  # magnitude matrix (using individual posterior modes here)
  mat <- log10(exp(jp$modes)+1);
  # weighted distance
  mode.fail.dist <- as.dist(1-do.call(rbind,mclapply(cell.names,function(nam1) {
    unlist(lapply(cell.names,function(nam2) {
      corr(cbind(mat[,nam1],mat[,nam2]),w=sqrt(sqrt(matw[,nam1]*matw[,nam2])))
    }))
  },mc.cores=n.cores)),upper=F);
  #    rownames(mode.fail.dist) <- colnames(gene.counts)
  #    colnames(mode.fail.dist) <- colnames(gene.counts)
  mode.fail.dist
}

my.plot.callback <- function(x) {
  library(rgl)
  plot3d(x, size=10, col=cols)
}

run.tsne <- function(my.dist, plot.callback=my.plot.callback, k=3, max_iter=5000, perplexity=5,initial_dims=10,...) {
  library(tsne)
  my.tsne <- tsne(my.dist, k=k, epoch_callback=plot.callback, max_iter=max_iter, perplexity = perplexity,whiten = T,initial_dims = initial_dims)
  rownames(my.tsne) <- labels(my.dist)
  my.tsne <- data.table(my.tsne,keep.rownames = T)
  my.tsne
}

do.dbscan <- function(data,eps){
  # assumes data is a data.table with column 1=row identifier
  # and column 2:end = data points
  my.dbscan <- copy(data)
  library(fpc)
  dbs=dbscan(my.dbscan[,2:ncol(my.dbscan),with=F],eps=eps,method="raw")
  my.dbscan[,cluster:=dbs$cluster]
  my.dbscan
}

colorby.gene <- function(data,xvar,yvar,colorvar){
  xcol=which(colnames(data)==xvar)
  ycol=which(colnames(data)==yvar)
  colorcol=which(colnames(data) %in% colorvar)
  cols=vecs2rgb(data[,colorcol])
  plt=plot(data[,c(xcol,ycol)],col=cols)
  plt
}
do.wss <- function(data,kmax=20,iter=20){
  wss=rep(0,20)
  for(i in 1:20){
    km <- kmeans(data,i,iter.max = 200)
    wss[i]=sum(km$withinss)
  }
  plt=plot(log10(wss),main="Within-group sum-of-squares",ylab=expression(paste(log[10]," wss")),xlab="k")
  plt=plt+lines(log10(wss))
  plt
}


plot.sequential <- function(geneset, counts.log, my.tsne, cols=NULL) {
  for(gene in geneset) {
    invisible(readline(prompt="Press [enter] to continue"))
    if(!is.na(match(gene, rownames(counts.log)))) {
      #        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
      #        panc.genes.max <- panc.genes.count
      panc.g.col <- counts.log[gene,]

      #        sum(panc.genes.max <= 1)
      #        panc.genes.max[panc.genes.max <= 1] <- 1
      #        panc.genes.max[panc.genes.max == 0] <- 0.0001
      #        panc.g.col <- as.numeric(log2(panc.genes.max))
      panc.g.col <- (panc.g.col-min(panc.g.col))/(max(panc.g.col)-min(panc.g.col))
      if(is.null(cols)) {
        cols <- colorRampPalette(c("Black", "Blue", "Yellow"))(100)
      }
      panc.g.col <- cols[panc.g.col*(length(cols)-2)+1]
      #        plot3d(my.rtsne$Y[,1], my.rtsne$Y[,2], my.rtsne$Y[,3], col=panc.g.col, size=10, main=gene)
      plot3d(my.tsne, col=panc.g.col, size=10, main=gene)
    }
    else {
      cat(gene, " not found\n");
    }
  }
}

geneset.colors <- function(gene, counts.log) {
  if(!is.na(match(gene[1], rownames(counts.log)))) {
    #        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
    #        panc.genes.max <- panc.genes.count
    panc.g.col1 <- counts.log[gene[1],]
    panc.g.col1 <- (panc.g.col1-min(panc.g.col1))/(max(panc.g.col1)-min(panc.g.col1))
    panc.g.col2 <- rep(0, dim(counts.log)[2])
    panc.g.col3 <- rep(0, dim(counts.log)[2])
    if(length(gene) > 1) {
      cat("two colors!")
      panc.g.col2 <- counts.log[gene[2],]
      panc.g.col2 <- (panc.g.col2-min(panc.g.col2))/(max(panc.g.col2)-min(panc.g.col2))
    }
    if(length(gene) > 2) {
      cat("three colors!")
      panc.g.col3 <- counts.log[gene[3],]
      panc.g.col3 <- (panc.g.col3-min(panc.g.col3))/(max(panc.g.col3)-min(panc.g.col3))
    }
    #        sum(panc.genes.max <= 1)
    #        panc.genes.max[panc.genes.max <= 1] <- 1
    #        panc.genes.max[panc.genes.max == 0] <- 0.0001
    #        panc.g.col <- as.numeric(log2(panc.genes.max))
    col <- rgb(panc.g.col1, panc.g.col2, panc.g.col3)
    #        cat(col)
    return(col)
  }
  return(NA)
}

vecs2rgb <- function(r, g=NULL, b=NULL) {
  r <- (r-min(r))/(max(r)-min(r))
  G <- rep(0, length(r))
  B <- rep(0, length(r))
  if(!is.null(g)) {
    G <- (g-min(g))/(max(g)-min(g))
  }
  if(!is.null(b)) {
    B <- (b-min(b))/(max(b)-min(b))
  }
  rgb(r, G, B)
}

darkvecs2rbg <- function(r, g=NULL, b=NULL,offset=0.2) {
  r <- (r-min(r))/(max(r)-min(r))
  G <- rep(0, length(r))
  B <- rep(0, length(r))
  if(!is.null(g)) {
    G <- (g-min(g))/(max(g)-min(g))
  }
  if(!is.null(b)) {
    B <- (b-min(b))/(max(b)-min(b))
  }
  rgb(r+offset*(1-r), G+offset*(1-G), B+offset*(1-B))
}

plot.coreg <- function(geneset, counts.log, my.tsne) {
  library(rgl)
  for(gene in geneset) {
    #    invisible(readline(prompt="Press [enter] to continue"))
    if(!is.na(match(gene[1], rownames(counts.log))) & !is.na(match(gene[2], rownames(counts.log)))) {
      #        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
      #        panc.genes.max <- panc.genes.count
      panc.g.col1 <- counts.log[gene[1],]
      panc.g.col1 <- (panc.g.col1-min(panc.g.col1))/(max(panc.g.col1)-min(panc.g.col1))
      panc.g.col2 <- rep(0, dim(counts.log)[1])
      panc.g.col3 <- rep(0, dim(counts.log)[1])
      if(length(gene) > 1) {
        cat("two colors!")
        panc.g.col2 <- counts.log[gene[2],]
        panc.g.col2 <- (panc.g.col2-min(panc.g.col2))/(max(panc.g.col2)-min(panc.g.col2))
      }
      if(length(gene) > 2) {
        cat("three colors!")
        panc.g.col3 <- counts.log[gene[3],]
        panc.g.col3 <- (panc.g.col3-min(panc.g.col3))/(max(panc.g.col3)-min(panc.g.col3))
      }
      #        sum(panc.genes.max <= 1)
      #        panc.genes.max[panc.genes.max <= 1] <- 1
      #        panc.genes.max[panc.genes.max == 0] <- 0.0001
      #        panc.g.col <- as.numeric(log2(panc.genes.max))
      col <- rgb(panc.g.col1, panc.g.col2, panc.g.col3)
      #        plot3d(my.rtsne$Y[,1], my.rtsne$Y[,2], my.rtsne$Y[,3], col=panc.g.col, size=10, main=gene)
      plot3d(my.tsne, col=col, size=10, main=gene)
    }
    else {
      cat(gene, " not found\n");
    }
  }
}

make.html <- function(counts.log, geneset, my.tsne, file.suffix="_expression_in_CIRM.html") {
  for(gene in geneset) {
    #    invisible(readline(prompt="Press [enter] to continue"))
    if(!is.na(match(gene, rownames(counts.log)))) {
      #        panc.genes.count <- counts[gene,]/colSums(counts)*1000000
      #        panc.genes.max <- panc.genes.count
      panc.g.col <- counts.log[gene,]

      #        sum(panc.genes.max <= 1)
      #        panc.genes.max[panc.genes.max <= 1] <- 1
      #        panc.genes.max[panc.genes.max == 0] <- 0.0001
      #        panc.g.col <- as.numeric(log2(panc.genes.max))
      panc.g.col <- (panc.g.col-min(panc.g.col))/(max(panc.g.col)-min(panc.g.col))
      panc.g.col <- colorRampPalette(c("Black", "Blue", "Yellow"))(100)[panc.g.col*98+1]
      #        plot3d(my.rtsne$Y[,1], my.rtsne$Y[,2], my.rtsne$Y[,3], col=panc.g.col, size=10, main=gene)
      #        plot3d(my.tsne, col=panc.g.col, size=10, main=gene)
      make.3dplot(my.tsne, panc.g.col, paste(gene, file.suffix, sep=""))

    }
    else {
      cat(gene, " not found\n");
    }
  }
}

make.duplicate.predictor <- function(counts) {
  library(randomForest)
  n.samples <- dim(counts)[2]
  cat("Number of samples: ", n.samples)
  n.genes <- dim(counts)[1]
  cat("Number of genes: ", n.genes)
  mixed.df <- as.data.frame(lapply(1:n.samples, function(x) {
    cat("x: ", x)
    r <- sample(1:n.samples, 2, replace=FALSE)
    merge.cells.wercc(counts[,r[1]], counts[,r[2]])
  }))
  both.df <- cbind(counts, mixed.df)
  both.df.norm <- norm.log.counts(both.df)
  rownames(both.df.norm) <- rownames(counts)
  y <- as.factor(c(rep(1,n.samples ), rep(2, n.samples)))
  my.rf1 <- randomForest(x=t(both.df.norm), y=y, ntree=2000, importance=TRUE, do.trace=T)
  cat(summary(my.rf1))
  my.rf1
}

merge.cells <- function(x, y) {
  # pick equal number of random counts from x and y
  names <- 1:length(x)
  ox <- x[names] > 0
  lx <- unlist(lapply(names[ox], function(i) {rep(i, x[i])}))

  oy <- y[names] > 0
  ly <- unlist(lapply(names[oy], function(i) {rep(i, y[i])}))

  num.to.pick <- mean(length(lx), length(ly))

  lnew <- c(sample(lx, as.integer(num.to.pick/2)), sample(lx, as.integer(num.to.pick/2)))
  zt <- table(lnew)
  z <- x-x
  z[as.integer(names(zt))] <- as.integer(zt)
  z
}

merge.cells.wercc <- function(xorig, yorig) {
  # pick equal number of random counts from x and y
  ercc.o <- grepl("^ERCC", names(xorig))
  x.ercc <- xorig[ercc.o]
  z <- xorig-xorig
  z[!ercc.o] <- xorig[!ercc.o]+yorig[!ercc.o]
  z[ercc.o] <- xorig[ercc.o]
  z
}

pick.clusters <- function(layout, classification) {
  plot(layout, col=classification, pch=19, cex=0.8)
  #    invisible(readline(prompt="Press [enter] to pick clusters"))
  selected.clusts <- unique(classification[identify(layout, labels=classification)])
  #    invisible(readline(prompt="Press [enter] to show selected cells"))
  selected.cells <- !is.na(match(classification, selected.clusts))
  plot(layout, col=selected.cells+1, pch=19, cex=0.8)
  selected.cells
}

# get distribution of counts and ERCC/total ratio
get.mm10.ercc.ratio <- function(data.counts){
  library(data.table)
  data.counts[,gene_type:="mm10"]
  data.counts[gene %like% "ERCC-",gene_type:="ERCC"]
  data.counts[,sum_counts:=sum(count),by=list(cell.name,gene_type)]

  cell.counts <- dcast.data.table(unique(data.counts[,list(cell.name,gene_type,sum_counts)]),cell.name ~ gene_type, value.var="sum_counts",fun.aggregate=mean,fill=0)
  cell.counts
}


get.gene.cor <- function(data, min.cells.expr = 1) {
  mat = t(data[, colSums(data > 0) > min.cells.expr])
  mat = mat - rowMeans(mat)
  cr = tcrossprod( mat/sqrt(rowSums(mat^2)) )

  return(cr)
}

filterByMinCor <- function(dataMatrix, n.genes = 2500, min.cells.expr = 2) {

  gene.min.cors <- get.min.cor(dataMatrix, min.cells.expr = min.cells.expr)
  genes.use <- names(sort(gene.min.cors)[1:n.genes])

  return(dataMatrix[, genes.use])
}

# return minimum (most negative) correlation value for each gene
get.min.cor <- function(data, min.cells.expr = 1) {
  library(matrixStats)

  cr = get.gene.cor(data, min.cells.expr)

  min.cor = colMins(cr, na.rm =T)
  names(min.cor) = colnames(cr)
  return(min.cor)
}

# return maximum correlation value for each gene
get.max.cor <- function(data, min.cells.expr = 1) {
  library(matrixStats)

  cr = get.gene.cor(data, min.cells.expr)
  diag(cr) <- NA # remove self correlations
  max.cor = colMaxs(cr, na.rm =T)
  names(max.cor) = colnames(cr)
  return(max.cor)
}

plot.cormap <- function(mat.cor,fname,width=12,height=12,method="complete"){
  nrows=nrow(mat.cor)
  hc.cor <- hclust(dist(mat.cor), method=method)
  hr.cor <- hclust(dist(mat.cor), method=method)
  cor.palette <- colorRampPalette(c("dodgerblue4", "lightsteelblue3", "white", "lightpink3", "red4"), space="Lab")
  cor.palette.breaks <- seq(-1,1,0.01)

  source("~/singlecell-pipeline/heatmap.ticks.R")
  pdf(fname, width=width, height=height)
  heatmap.ticks(mat.cor,trace="none",density="none",
                col=cor.palette, breaks=cor.palette.breaks,
                Colv=as.dendrogram(hc.cor), Rowv=as.dendrogram(hr.cor),
                dendrogram="column",labCol=F, key.title=NA, keysize=1,
                key.xlab="Pearson correlation", cexRow=1.5*sqrt(height/nrows),mar=c(5,8))
  dev.off()
}

# submsample number of genes detected per cell
# do not subsample cells with fewer genes than minimum
subsample.numgenes <- function(data, ng.min = 4500, n.cores = 8, from.mat = T, max.frac.expr = 0.95,
                               normalize = T){

  # munging
  if(!from.mat){
    mat.data = cast.to.df(data, to.matrix = T)
  } else {
    mat.data = data
  }

  CELLS <- rownames(mat.data)
  genes <- colnames(mat.data)

  print(summary(rowSums(mat.data > 0)))

  # list of genes to avoid subsampling because they are constitutively expressed
  genes.constitutive <- genes[colSums(mat.data > 0) > max.frac.expr * length(CELLS)]
  genes.use <- genes[!genes %in% genes.constitutive]
  if(length(genes.constitutive) > 0) mat.constitutive <- mat.data[, genes.constitutive]
  mat.use <- mat.data[, genes.use]
  registerDoMC(n.cores)


  set.rand.zeros <- function(x, ng.keep) {
    if(ng.keep < sum(x > 0)){
      y = numeric(length(x))
      rand.include = sample(1:length(y), ng.keep)
      y[rand.include] <- x[rand.include]
      return(y)
    }
    else{
      return(x)
    }

  }

  # subsample
  data.subsamp = foreach(cell = rownames(mat.use), .combine = rbind) %dopar% {

    cell.expr  = mat.use[cell, ]
    g.expr     = names(cell.expr[cell.expr > 0])
    ng.current = sum(cell.expr > 0)

    # possibility of non-hard subsample: weighted average of
    # current and target number of genes
    # ng.subsamp = ceiling(0*length(g.expr) + 1*ng.min)

    if(ng.current > ng.min){
      g.remove = sample(g.expr, ng.current - ng.min)
      cell.expr[g.remove] <- 0
    }

    cell.expr
  }
  rownames(data.subsamp) <- rownames(mat.use)

  if(length(genes.constitutive) > 0) data.subsamp <- cbind(data.subsamp, mat.constitutive)

  if(normalize) data.subsamp <- log10( 1 + 1e6*data.subsamp/rowSums(data.subsamp) )

  print(summary(rowSums(data.subsamp > 0)))

  if(!from.mat){
    data.out = data.table(melt(data.subsamp, varnames = c('cell.name', 'gene'), value.name = 'log10.cpm'))
    return(data.out)
  } else {
    return(data.subsamp)
  }


}


# data is a matrix with rows = cells and columns = genes, log10.cpm
# cell.info is a data.table
dimplots.pca.lowMem <- function(data,cell.info,dir,suffix="",ncp=10,plot.sigs=T,max.genes=30,min.cells.detect=1,cexCol=0.6,scramble=F,expression.values="log10.cpm",genes.sig=NA,cex.sig=1,
                                plot.gene.cor=F,sig.exptCol=T,n.cores=detectCores(),plot.sigHeatmaps=F,PCA.type = 'robust',min.numgenes = NA,
                                ng.subsample = NA,...){
  if(!file.exists(dir)) system(paste0('mkdir ', dir))
  if(!file.exists(file.path(dir, 'biplots/'))) system(paste0('mkdir ', file.path(dir, 'biplots/')))

  if(is.na(genes.sig)){
    genes.sig=max.genes
  }

  # make sure directory name ends in slash (need to change to run on Windows!)
  if(substr(dir,nchar(dir),nchar(dir))!="/"){
    dir <- paste0(dir,"/")
  }
  # make named character to add custom experiment colors to ggplots
  setkey(cell.info)
  expt.col <- unique(cell.info[,.(experiment,experiment_color)]); setkey(expt.col, experiment)
  expt2col <- expt.col[['experiment_color']]; names(expt2col) <- expt.col[['experiment']]
  expt_cols=data.frame(row.names = cell.info[,cell.name],experiment_color=cell.info[,experiment_color])


  # remove cells with less than minimum number of genes (if specified)
  if(!is.na(min.numgenes)){
    data[, num.genes := sum(get(expression.values) > 0), by = 'cell.name']
    data = data[num.genes > min.numgenes]
  }

  # only subsample the number of genes for calculating PCs.
  # Keep all genes for calculating scores
  if(!is.na(ng.subsample)){
    data.pca = subsample.numgenes(data, ng.min = ng.subsample, from.mat = F)
  }
  else {
    data.pca = data
  }

  print("calculating PCA")
  mat.pca = cast.to.df(data.pca, expression.values = expression.values, to.matrix = T)
  if(PCA.type == 'robust'){
    if(scramble){pcalist <- do.robpca.scramble(mat.pca,ncp = ncp,expression.values=expression.values); suffix="_scrambled"}
    else{pcalist <- do.robpca(mat.pca,ncp = ncp,expression.values=expression.values,...)}

  } else if(PCA.type == 'classical'){
    pcalist <- do.pca(mat.pca,ncp = ncp,expression.values=expression.values,...)
  }

  loadings <- pcalist[[3]]
  cast.pca <- pcalist[[2]]
  rm(pcalist);gc()
  setnames(cast.pca, c(paste0('PC',1:ncp),'cell.name'))
  setnames(loadings, c(paste0('PC',1:ncp),'gene'))

  # write loadings and scores to a csv
  # write.csv(cast.pca,file=paste0(dir,"PCscores",suffix,".csv"),quote=F,row.names=T)
  write.csv(loadings,file=paste0(dir,"PCloadings",suffix,".csv"),row.names=T)
  gc()
  cell.info <- merge(cell.info, cast.pca, by="cell.name")
  #data.pca[,numcells_det:=sum(log10.cpm>0),by=gene]; data.pca <- data.pca[numcells_det > min.cells.detect]
  dims <- colnames(cast.pca)[1:(ncol(cast.pca)-1)]


  out.sig <- data.table(cell.name=cell.info[,cell.name])
  out.sigscore <- data.table(cell.name=cell.info[,cell.name])

  gc()
  # calculate PC signatures as the sum of the top 60 genes by PC
  print("calculating modified PC scores and saving plots")
  library(foreach)
  pc.allscores <- foreach(dim=dims,.combine='cbind') %dopar% {

    setorderv(loadings,dim,-1); genes.pos <- loadings[,gene][1:max.genes]
    setorderv(loadings,dim,1); genes.neg <- loadings[,gene][max.genes:1]

    genes.use <- unique(c(genes.pos,genes.neg))
    data.genes.use <- data[gene %in% genes.use]

    # following code needed to prevent cells with no expression of genes.pos or genes.neg being removed from the dataset
    setkey(data)
    cells.noExpr=cell.info[!cell.name %in% unique(data.genes.use[,cell.name])]
    if(nrow(cells.noExpr) > 0){
      cells.noExpr[,gene:=genes.pos[1]];cells.noExpr[,log10.cpm:=0];
      data.genes.use=rbindlist(list(data.genes.use[,.(cell.name,gene,log10.cpm)],cells.noExpr),use.names = T,fill = F)
      setkey(data)
      print(length(unique(data.genes.use[gene %in% genes.use,cell.name])))

    }

    # generate signature value of top PC +/- correlated genes
    cast.sig <- get.sig(data.genes.use, genes.pos,genes.neg, from.mat = F)
    cast.sig <- merge(cast.sig, unique(cell.info[,.(cell.name,experiment,experiment_color)]),by='cell.name')


    # plot signature biplots and histograms
    if(sig.exptCol){
      plt=ggplot(cast.sig,aes(x=pos,y=neg,color=experiment))+geom_point()+
        scale_color_manual(values=expt2col)+theme_bw()+coord_fixed()
      ggsave(plt,file=paste0(dir,dim,"_sigScatter",suffix,".pdf"),height=4,width=8)
    }
    else if(!sig.exptCol){
      plt=ggplot(cast.sig,aes(x=pos,y=neg))+geom_point()+
        theme_bw()+coord_fixed()
      ggsave(plt,file=paste0(dir,dim,"_sigScatter",suffix,".pdf"),height=4,width=6)

    }


    # plot heatmap of top 15 signature genes ordered by signature value
    if(plot.sigHeatmaps){
      setorderv(loadings,dim,-1); genes.pos.sig <- loadings[,gene][1:15]
      setorderv(loadings,dim,1); genes.neg.sig <- loadings[,gene][15:1]
      plot.sigHeatmap.wExptCols(data,genes.pos,genes.neg,fname = paste0(dir,dim,"_sigHeatmap",suffix,".pdf"),width=15,height=10,cexCol=cexCol)
      # same heatmap with genes ordered by hierarchical clustering
      cast.sig[,frac.pos:=pos-neg]; setorder(cast.sig,frac.pos)
      genes.use <- c(genes.pos.sig,genes.neg.sig)
      cast.plot <- cast.to.df(data[gene %in% genes.use],annot='experiment_color');
      cast.plot <- cast.plot[match(cast.sig[['cell.name']],rownames(cast.plot)),c('experiment_color',genes.use)]
      cellorder.heatmap.wExptCols(cast.plot,fname = paste0(dir,dim,"_geneCluster",suffix,".pdf"),width=15,height=10,cexCol=cex.sig)

      # genes and cells ordered by HC
      cell.heatmap.wExptCols(cast.plot, fname = paste0(dir,dim,"_","hclust",suffix,".pdf"),width=15,height=10,cexCol=cex.sig)

    }

    gc()

    # plot 2d signature density
    # plt=ggplot(cast.sig,aes(x=pos,y=neg))+stat_bin2d(bins=10)+scale_fill_gradientn(colours=scale.farmhouse.rose)
    # ggsave(plt,file=paste0(dir,dim,"_sigDensity",suffix,".pdf"),height=4,width=5)
    #
    # plot signature histogram
    plt=ggplot(cast.sig,aes(pos-neg,..density..))+geom_histogram(binwidth=0.05,fill='grey50')+ geom_density(color='white')+theme_bw()
    ggsave(plt,file=paste0(dir,dim,"_sigHist",suffix,".pdf"),height=4,width=5)

    gc()
    cast.sig[,pc.score:=pos-neg];
    # cast.sig[,pc.score.mc:=pos.mc-neg.mc]
    cast.out.sig <- cast.sig[,.(cell.name,pos,neg,pc.score)]#,pos.mc,neg.mc,pc.score.mc)]
    setnames(cast.out.sig,c('cell.name',paste0(dim,'.pos'),paste0(dim,'.neg'),paste0(dim,'.score')))
    # paste0(dim,'.pos.mc'),paste0(dim,'.neg.mc'),paste0(dim,'.score.mc')))
    setorder(cast.out.sig,cell.name);cast.out.sig[,cell.name:=NULL]

    as.matrix(cast.out.sig)

  }

  pc.allscores=data.table(pc.allscores,keep.rownames = F)
  setorder(cell.info,cell.name)
  pc.allscores[,cell.name:=cell.info[,cell.name]]
  # print(pc.allscores)
  pc.allscores = merge(pc.allscores,cast.pca,by='cell.name')
  write.csv(pc.allscores,file=paste0(dir,"PC_allscores",suffix,".csv"),quote=F,row.names=T)

  list(cast.pca)
}


barb.cormap <- function(mat.cor,fname,width=12,height=12,method="complete",cex=0.5,mincor=-1,maxcor=1){
  library(lattice)
  library(cba)
  rowdist <- dist(mat.cor)
  coldist <- dist(mat.cor, by_rows = F)
  hc.cor <- hclust(coldist, method=method)
  hr.cor <- hclust(rowdist, method=method)

  optimal.row <- order.optimal(rowdist,hr.cor$merge)
  optimal.col <- order.optimal(coldist,hc.cor$merge)

  ord.row <- optimal.row$order
  ord.col <- optimal.col$order

  plt = levelplot(mat.cor[ord.row,ord.col],xlab=NULL,ylab=NULL,
                  at=do.breaks(c(mincor-0.01,maxcor+0.01),19),scales=list(x=list(rot=90),cex=cex),
                  colorkey=list(space="top"),
                  col.regions=colorRampPalette(c("dodgerblue4", "dodgerblue", "white", "lightcoral", "firebrick4"), space="Lab"))
  pdf(fname,width=width,height=height)
    print(plt)
  dev.off()

  return(list(rownames(mat.cor[ord.row, ]), colnames(mat.cor[, ord.col]), plt, hc.cor, hr.cor))
}

cell.heatmap <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,method="ward",cexRow=1,cexCol=1){
  data <- data[which(apply(data,1,sd)>0),]
  data <- data[,which(apply(data,2,sd)>0)]
  source("~/singlecell-pipeline/heatmap.ticks.R")
  coldist=as.dist(1-cor(data,method="spearman"));
  rowdist=as.dist(1-cor(t(data),method="spearman"))
  hc <- hclust(coldist, method=method)
  hr <- hclust(rowdist, method=method)

  library(cba)
  optimal.row <- order.optimal(rowdist,hr$merge)
  optimal.col <- order.optimal(coldist,hc$merge)
  hr$merge <- optimal.row$merge; hr$order <- optimal.row$order
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order

  #color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  color.palette = colorRampPalette(c("blue4","white","red3"), space="Lab")
  palette.breaks <- seq(min(data),max(data), (max(data)-min(data))/21)

  pdf(fname,width=width,height=height)
  # par(bg='grey70')
  heatmap.ticks(as.matrix(data),trace="none",density="none",
                Colv = as.dendrogram(hc),Rowv = as.dendrogram(hr),
                key=T,scale="none",dendrogram="both",col=color.palette,breaks=palette.breaks,
                labRow=rownames(data), labCol=colnames(data),
                cexRow=cexRow, cexCol=cexCol, keysize=0.8, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(round(min(palette.breaks),digits = 1),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                mar=c(12,18))
  dev.off()

  return(list(hc,hr))
}

cell.heatmap.wExptCols <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,method="ward",cexRow=1,cexCol=1){
  library(cba)
  data.noex <- data[,2:ncol(data)]
  data.noex <- data.noex[which(apply(data.noex,1,sd)>0),]
  data.noex <- data.noex[,which(apply(data.noex,2,sd)>0)]
  data <- data[rownames(data.noex),c("experiment_color",colnames(data.noex))]
  source("~/singlecell-pipeline/heatmap.ticks.R")
  coldist=as.dist(1-cor(data[,2:ncol(data)],method="pearson"));
  rowdist=as.dist(1-cor(t(data[,2:ncol(data)]),method="pearson"))
  hc <- hclust(coldist, method=method)
  hr <- hclust(rowdist, method=method)

  optimal.row <- order.optimal(rowdist,hr$merge)
  optimal.col <- order.optimal(coldist,hc$merge)
  hr$merge <- optimal.row$merge; hr$order <- optimal.row$order
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order

  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  palette.breaks <- seq(min(data[,2:ncol(data)]),max(data[,2:ncol(data)]), (max(data[,2:ncol(data)])-min(data[,2:ncol(data)]))/21)

  pdf(fname,width=width,height=height)
  par(bg='grey80')
  heatmap.ticks(as.matrix(data[,2:ncol(data)]),trace="none",density="none",
                Colv = as.dendrogram(hc),Rowv = as.dendrogram(hr),
                key=T,scale="none",dendrogram="both",col=color.palette,breaks=palette.breaks,
                labRow=NA, labCol=colnames(data[,2:ncol(data)]),
                cexRow=cexRow, cexCol=cexCol, keysize=0.8, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                mar=c(12,10),RowSideColors=as.character(data[,1]))
  dev.off()

  return(list(hc,hr))
}

# for use with expression normalized by gene
norm.heatmap.wExptCols <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,method="ward",cexRow=1,cexCol=1){
  data.noex <- data[,2:ncol(data)]
  data.noex <- data.noex[which(apply(data.noex,1,sd)>0),]
  #   print(data.noex[1:10,1:10])
  data.noex <- data.noex[,which(apply(data.noex,2,sd)>0)]
  data <- data[rownames(data.noex),c("experiment_color",colnames(data.noex))]
  source("~/singlecell-pipeline/heatmap.ticks.R")
  coldist=as.dist(1-cor(data[,2:ncol(data)],method="spearman"));
  rowdist=as.dist(1-cor(t(data[,2:ncol(data)]),method="spearman"))
  hc <- hclust(coldist, method=method)
  hr <- hclust(rowdist, method=method)

  library(cba)
  optimal.row <- order.optimal(rowdist,hr$merge)
  optimal.col <- order.optimal(coldist,hc$merge)
  hr$merge <- optimal.row$merge; hr$order <- optimal.row$order
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order
  #
  color.palette = colorRampPalette(c("darkgreen","white","firebrick"), space="Lab")

  mode.data <- as.numeric(names(tail(sort(table(unlist(data[,2:ncol(data)]))),1)))
  max.data <- max(data[,2:ncol(data)]); min.data <- min(data[,2:ncol(data)]);
  palette.breaks <- c(seq(min.data,mode.data,(mode.data-min.data)/9),mode.data,seq(mode.data,max.data,(max.data-mode.data)/10))

  pdf(fname,width=width,height=height)
  heatmap.ticks(as.matrix(data[,2:ncol(data)]),trace="none",
                Colv = as.dendrogram(hc),Rowv = as.dendrogram(hr),
                key=T,scale="none",dendrogram="both",col=color.palette,breaks=palette.breaks,
                labRow=NA, labCol=colnames(data[,2:ncol(data)]),
                cexRow=cexRow, cexCol=cexCol, keysize=0.8, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                mar=c(12,10),RowSideColors=as.character(data[,1]))
  dev.off()

  return(list(hc,hr))
}

cell.heatmap.wGeneCols <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,method="ward",cexRow=1,cexCol=1){
  data.noex <- data[,2:ncol(data)]
  data.noex <- data.noex[,which(apply(data.noex,2,sd)>0)]
  data <- data[,c("gene_color",colnames(data.noex))]
  source("~/singlecell-pipeline/heatmap.ticks.R")
  coldist=as.dist(1-cor(data[,2:ncol(data)],method="spearman"));
  rowdist=as.dist(1-cor(t(data[,2:ncol(data)]),method="spearman"))
  hc <- hclust(coldist, method=method)
  hr <- hclust(rowdist, method=method)

  library(cba)
  optimal.row <- order.optimal(rowdist,hr$merge)
  optimal.col <- order.optimal(coldist,hc$merge)
  hr$merge <- optimal.row$merge; hr$order <- optimal.row$order
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order
  #
  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  palette.breaks <- seq(min(data[,2:ncol(data)]),max(data[,2:ncol(data)]), (max(data[,2:ncol(data)])-min(data[,2:ncol(data)]))/21)

  pdf(fname,width=width,height=height)
  heatmap.ticks(as.matrix(data[,2:ncol(data)]),trace="none",density="none",
                Colv = as.dendrogram(hc),Rowv = as.dendrogram(hr),
                key=T,scale="none",dendrogram="both",col=color.palette,breaks=palette.breaks,
                labRow=NA, labCol=colnames(data[,2:ncol(data)]),
                cexRow=cexRow, cexCol=cexCol, keysize=0.8, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                mar=c(12,10),RowSideColors=as.character(data[,1]))
  dev.off()

  return(list(hc,hr))
}

# order cells by input order and genes (columns) by clustering
cellorder.heatmap <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,method="ward",cexRow=1,cexCol=1){
  data <- data[,which(apply(data,2,sd)>0)]
  source("~/singlecell-pipeline/heatmap.ticks.R")
  coldist=as.dist(1-cor(data,method='spearman'));
  hc <- hclust(coldist, method=method)

  library(cba)
  optimal.col <- order.optimal(coldist,hc$merge)
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order

  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  #color.palette = colorRampPalette(c('red','red3','black','green3','green'), space="Lab")
  palette.breaks <- seq(min(data),max(data), (max(data)-min(data))/41)

  pdf(fname,width=width,height=height)
  heatmap.ticks(as.matrix(data),trace="none",density="none",
                Colv = as.dendrogram(hc),Rowv = F,
                key=T,scale="none",dendrogram="column",col=color.palette,breaks=palette.breaks,
                labRow=rownames(data), labCol=colnames(data),
                cexRow=cexRow, cexCol=cexCol, keysize=1.2, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                mar=c(12,15))
  dev.off()
}

heatmap.pal <- c("#18176D","#723DC9","white","#B7D71F","#65BA00")
# first column of data must be the color with column title "experiment_color"
cellorder.heatmap.wExptCols <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,method="ward",cexRow=1,cexCol=1){
  data.noex <- data[,2:ncol(data)]
  #   data.noex <- data.noex[which(apply(data.noex,1,sd)>0),]
  data.noex <- data.noex[,which(apply(data.noex,2,sd)>0)]
  #   data <- data[rownames(data.noex),c("experiment_color",colnames(data.noex))]
  data <- data[,c("experiment_color",colnames(data.noex))]
  source("~/singlecell-pipeline/heatmap.ticks.R")
  coldist=as.dist(1-cor(data[,2:ncol(data)],method='spearman'));
  hc <- hclust(coldist, method=method)

  library(cba)
  optimal.col <- order.optimal(coldist,hc$merge)
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order

  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  # color.palette = colorRampPalette(heatmap.pal, space="Lab")
  #   col.gen=colorRampPalette( brewer.pal(11,'PRGn')[11:1],space="Lab")
  #   color.palette = colorRampPalette( c(col.gen(100)[91],col.gen(100)[82],col.gen(100)[50],
  #                                       col.gen(100)[18],col.gen(100)[10]), space="Lab")
  #   color.palette = colorRampPalette(c('black','red3','white'), space="Lab")
  palette.breaks <- seq(min(data.noex),max(data.noex), (max(data.noex)-min(data.noex))/21)

  pdf(fname,width=width,height=height)
  heatmap.ticks(as.matrix(data[,2:ncol(data)]),trace="none",density="none",
                Colv = as.dendrogram(hc),Rowv = F,
                key=T,scale="none",dendrogram="column",col=color.palette,breaks=palette.breaks,
                labRow=NA, labCol=colnames(data[,2:ncol(data)]),
                cexRow=cexRow, cexCol=cexCol, keysize=1.2, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                RowSideColors=as.character(data[,1]),
                mar=c(12,15))
  dev.off()
}

# data should be a data.frame (not data.table)
no.dendro.heatmap <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,method="ward",cexRow=1,cexCol=1){
  source("~/singlecell-pipeline/heatmap.ticks.R")

  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  palette.breaks <- seq(min(data),max(data), (max(data)-min(data))/21)
  pdf(fname,width=width,height=height)
  heatmap.ticks(as.matrix(data),trace="none",density="none",
                Colv = F,Rowv = F,
                key=T,scale="none",dendrogram="none",col=color.palette,breaks=palette.breaks,
                labRow=NA, labCol=colnames(data),
                cexRow=cexRow, cexCol=cexCol, keysize=1.2, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=0.8, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                mar=c(12,10))
  dev.off()
}

no.dendro.heatmap.wExptCols <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,cexRow=1,cexCol=1){
  source("~/singlecell-pipeline/heatmap.ticks.R")

  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  palette.breaks <- seq(min(data[,2:ncol(data)]),max(data[,2:ncol(data)]), (max(data[,2:ncol(data)])-min(data[,2:ncol(data)]))/21)
  pdf(fname,width=width,height=height)

  heatmap.ticks(as.matrix(data[,2:ncol(data)]),trace="none",density="none",
                Colv = F,Rowv = F,
                key=T,scale="none",dendrogram="none",col=color.palette,breaks=palette.breaks,
                labRow=NA, labCol=colnames(data[,2:ncol(data)]),
                cexRow=cexRow, cexCol=cexCol, keysize=1.2, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                RowSideColors=as.character(data[,1]),
                mar=c(12,10))
  dev.off()
}

# plot heatmap of data, rows (cells) ordered by expression signature of two opposing gene sets and
# columns (genes) ordered by input vectors
plot.sigHeatmap.wExptCols <- function(data,genes.pos,genes.neg,fname,cluster=F,cexCol=0.6,controlPlot=F,genes.ctl=NULL,...){
  cast.sig <- get.sig(data,genes.pos=genes.pos,genes.neg=genes.neg,write_csv = F)
  cast.sig[,sig:=pos-neg]; setorder(cast.sig,sig)

  if(controlPlot){
    genes.use <- genes.ctl
    cast.plot <- cast.to.df(data[gene %in% genes.use],annot='experiment_color');
    cast.plot <- cast.plot[match(cast.sig[['cell.name']],rownames(cast.plot)),c('experiment_color',genes.use)]


    if(cluster){cellorder.heatmap.wExptCols(cast.plot,fname,cexCol=cexCol,...)}
    else{no.dendro.heatmap.wExptCols(cast.plot,fname,...)}
  }
  else{
    genes.use <- c(genes.pos,genes.neg)
    cast.plot <- cast.to.df(data[gene %in% genes.use],annot='experiment_color');
    cast.plot <- cast.plot[match(cast.sig[['cell.name']],rownames(cast.plot)),c('experiment_color',genes.pos,genes.neg)]


    if(cluster){cellorder.heatmap.wExptCols(cast.plot,fname,cexCol=cexCol,...)}
    else{no.dendro.heatmap.wExptCols(cast.plot,fname,...)}
  }

  # cast.plot

}



sigorder.heatmap.wExptCols <- function(data,fname,width=width,height=height,rowCols=NA,colCols=NA,cexRow=1,cexCol=1){
  source("~/singlecell-pipeline/heatmap.ticks.R")

  color.palette = colorRampPalette(c("midnightblue","dodgerblue3","white","goldenrod1","darkorange2"), space="Lab")
  palette.breaks <- seq(min(data[,2:ncol(data)]),max(data[,2:ncol(data)]), (max(data[,2:ncol(data)])-min(data[,2:ncol(data)]))/21)
  #   print(palette.breaks)
  #   print(data)
  pdf(fname,width=width,height=height)
  heatmap.ticks(as.matrix(data[,2:ncol(data)]),trace="none",density="none",
                Colv = F,Rowv = F,
                key=T,scale="none",dendrogram="none",col=color.palette,breaks=palette.breaks,
                labRow=NA, labCol=colnames(data[,2:ncol(data)]),
                cexRow=cexRow, cexCol=cexCol, keysize=1.2, key.title=NA, key.xlab=expression(paste(log[10], " cpm")),
                key.xlab.size=2, key.ticks = c(min(palette.breaks),round(mean(palette.breaks),digits = 1),round(max(palette.breaks),digits = 1)),
                RowSideColors=as.character(data[,1]),
                mar=c(12,10))
  dev.off()
}
do.pca <- function(data.counts,min_cells=0,min_sd=0,ncp=10,expression.values="log10.cpm", from.mat = T){
  library(FactoMineR)
  if(!from.mat){
    if(min_cells>0){
      data<-filter.genes(data.counts,expression.values,min_sd,min_cells)

    } else{ cast.counts <- cast.to.df(data.counts,expression.values) }
  } else {
    cast.counts = data.counts[, colSums(data.counts > 0) > min_cells]
  }

  res <- PCA(cast.counts, ncp=ncp, graph=F,scale.unit = F)
  coord <- as.data.frame(res$ind$coord[, 1:ncp])
  coord$cell.name <- rownames(coord); coord <- data.table(coord)
  loadings <- as.data.frame(res$var$coord[, 1:ncp])
  loadings$gene <- rownames(loadings); loadings <- data.table(loadings)
  list(res,coord,loadings)
}

do.robpca <- function(data.counts,min_cells=0,min_sd=0,ncp=10,expression.values="log10.cpm",from.mat = T,...){
  library(rrcov)

  if(!from.mat){
    if(min_cells>0){
      data<-filter.genes(data.counts,expression.values,min_sd,min_cells)

    } else{ cast.counts <- cast.to.df(data.counts,expression.values) }
  } else {
      cast.counts = data.counts[, colSums(data.counts > 0) > min_cells]
  }

  #   data[,expr.mean:=mean(get(expression.value)),by=.(gene)]
  #   data[,scaled.exp:=get(expression.value)/expr.mean]

  print('starting robustPCA')
  pca <- PcaHubert(cast.counts, k = ncp,kmax=ncp,...)
  rm(cast.counts);gc()

  scores <- as.data.frame(getScores(pca));scores$cell.name <- rownames(scores); scores <- data.table(scores);setkey(scores)
  var   <- as.data.frame(getLoadings(pca)); var$gene <- rownames(var);var <- as.data.table(var);setkey(var)
  gc()
  list(pca, scores, var)
}
do.robpca.scramble <- function(data.counts,min_cells=0,min_sd=0,ncp=10,expression.values="log10.cpm",...){
  library(rrcov)
  cast.counts <- cast.to.df(data.counts,expression.values)
  for(j in 1:ncol(cast.counts)){
    cast.counts[,j] <- cast.counts[sample(nrow(cast.counts),replace = F),j]
  }
  if(min_cells>0){
    data<-filter.genes(data.counts,expression.values,min_sd,min_cells)

  }
  #   data[,expr.mean:=mean(get(expression.value)),by=.(gene)]
  #   data[,scaled.exp:=get(expression.value)/expr.mean]

  pca <- PcaHubert(cast.counts, k = ncp,kmax=ncp,...)

  scores <- as.data.frame(getScores(pca));scores$cell.name <- rownames(scores); scores <- data.table(scores);setkey(scores)
  var   <- as.data.frame(getLoadings(pca)); var$gene <- rownames(var);var <- as.data.table(var);setkey(var)
  gc()
  list(pca, scores, var)
}

# function to cast a melted data.table (columns cell.name, gene, log10.cpm)
# to a data.frame with rownames = cell.names and colnames = gene
cast.to.df <- function(data.counts,expression.values="log10.cpm",annot="none", to.matrix = F, genes.use = NULL){
  if(!is.null(genes.use)) data.counts = data.counts[gene %in% genes.use]

  if(annot != "none"){
    cast.counts <- dcast.data.table(data.counts, cell.name ~ gene, value.var=expression.values, fill = 0)
    setkey(data.counts)
    cast.counts <- merge(cast.counts,unique(data.counts[,.(cell.name,get(annot))]),by="cell.name")

    setnames(cast.counts,c(colnames(cast.counts)[1:(ncol(cast.counts)-1)],annot))
    cast.counts <- as.data.frame(cast.counts)
  }
  else{
    cast.counts <- as.data.frame(dcast.data.table(data.counts, cell.name ~ gene, value.var=expression.values,fill=0))
    gc()
  }

  rownames(cast.counts)<-cast.counts$cell.name; cast.counts <- cast.counts[,c(ncol(cast.counts),2:(ncol(cast.counts)-1))]
  if(to.matrix) cast.counts = as.matrix(cast.counts)
  cast.counts
}

# function to cast a melted data.table (columns cell.name, gene, log10.cpm)
# to a (numeric) matrix with rownames = cell.names and colnames = gene
# Matrices cannot be annotated as the purpose is for them to not library type-checking
# for speed and thus all elements must be numeric
cast.to.mat <- function(data.table, expression.values="log10.cpm", genes.use = NULL){


  cast.mat = matrix()
}

filter.genes <- function(data.counts,expression.value="log10.cpm",min_sd=0,min_cells=0){
  data=copy(data.counts);
  data[,numcells_exp:=sum(get(expression.value)>0),by=gene]; data<-data[numcells_exp>min_cells]
  data[,gene_sd:=sd(get(expression.value)),by=gene]; data<-data[gene_sd>min_sd]
  data
}

get.expt <- function(cell.name){
  expt= substr(cell.name,regexpr("BTN",cell.name),regexpr("BTN",cell.name)+4)
}

plot.insilico <- function(data.xy,data.counts,x,y,stainR,stainG=NULL,stainB=NULL,expression.values="log10.cpm",...){
  data.genes <- data.counts[gene %in% c(stainR,stainG,stainB),list(cell.name,gene,get(expression.values))]
  setnames(data.genes,c("cell.name","gene",expression.values))
  data.genes <- data.genes[cell.name %in% data.xy[,cell.name]]
  setkey(data.genes)
  cast.genes <- cast.to.df(data.genes, expression.values=expression.values)
  cast.genes <- cast.genes[sort(rownames(cast.genes)),]

  if (ncol(cast.genes)==1) {cols <- vecs2rgb(cast.genes[1])}
  else if (ncol(cast.genes)==2) {cols <- vecs2rgb(cast.genes[,stainR],cast.genes[,stainG])}
  else { cols <- vecs2rgb(cast.genes[,stainR],cast.genes[,stainG],cast.genes[,stainB])}
  names(cols) <- rownames(cast.genes)

  setkey(data.xy)
  coords <- as.data.frame(unique(data.xy[,list(cell.name,get(x),get(y))])); colnames(coords)<-c("cell.name",x,y)
  rownames(coords)<-coords$cell.name; coords <- coords[,2:3]
  cast.xy <- merge(cast.genes,coords,by="row.names"); rownames(cast.xy)<-cast.xy$Row.names; cast.xy<-cast.xy[,-1]
  cast.xy <- cast.xy[sort(rownames(cast.xy)),]
  cast.xy$cell.name <- rownames(cast.xy)
  colnames(cast.xy)<-c("x","y","cell.name")
  cast.xy <- cast.xy[,1:3]
  list(cast.xy, cols)
  #   plt=ggplot(cast.xy,aes(x=x,y=y,color=cell.name))+geom_point()+scale_color_manual(values=cols)+theme_bw()
  #   plt
}

dark.insilico <- function(cast.plot,x,y,red,green='',blue='',offset=0.3,fprefix='',height=5,width=5,psize=5,cex=2,axlim=T,...){
  cols = 'white'
  pdf(paste0(fprefix,'_',x,'_',y,'_',red,'_',green,'_',blue,'.pdf'),height=height,width=width)
  par(bg='black');

  if(axlim){
    plot(cast.plot[[x]],cast.plot[[y]],col=darkvecs2rbg(cast.plot[[red]],cast.plot[[green]],cast.plot[[blue]],offset=offset),cex=cex,
         xlab=x,ylab=y,col.lab=cols,bty='n',...)
  }
  else{
    plot(cast.plot[[x]],cast.plot[[y]],col=darkvecs2rbg(cast.plot[[red]],cast.plot[[green]],cast.plot[[blue]],offset=offset),cex=cex,
         xlab=x,ylab=y,col.lab=cols,bty='n',xlim=c(0,1),ylim=c(0,1),...)
  }
  axis(1,col=cols,col.axis=cols,col.ticks=cols);axis(2,col=cols,col.axis=cols,col.ticks=cols);box(col=cols,bty='l');

  dev.off()

  par(bg='black')
  print(plot(cast.plot[[x]],cast.plot[[y]],col=darkvecs2rbg(cast.plot[[red]],cast.plot[[green]],cast.plot[[blue]],offset=offset),cex=cex,
             xlab=x,ylab=y,col.lab=cols,bty='n',...))
  axis(1,col=cols,col.axis=cols,col.ticks=cols);axis(2,col=cols,col.axis=cols,col.ticks=cols);box(col=cols,bty='l');
}



make.mst <- function(dist.mat,vertex.color=NULL,vertex.size=5){
  library(igraph)
  adjc.mat<-copy(dist.mat)
  adjc.mat[lower.tri(adjc.mat, diag=TRUE)]<-0

  #Create the graph
  g_adjc <- graph.adjacency(adjc.mat > 0,weighted=TRUE,diag=FALSE,mode="upper")
  E(g_adjc)$weight<-round(t(adjc.mat)[(t(adjc.mat))>0],2)

  #Generate Minimum Spanning Tree
  mst_adjc <- minimum.spanning.tree(g_adjc)

  # Plot
  set.seed(12)
  plot(mst_adjc,vertex.size=vertex.size,vertex.color=vertex.color, vertex.label=NA)
  #hist(shortest.paths(mst_adjc), 1000, main="Histogram of shortest paths")
  #plot(density(shortest.paths(mst_adjc)), main="Density of shortest paths")

  list(g_adjc,mst_adjc)
}

group.mst <- function(mst){
  library(ggthemes, lib.loc = libloc)
  par(mfcol=c(4,4))
  #   for(i in c(10,20,30,40,50,60,70,90,92,94,96,98,100,102)){
  #     comm <- walktrap.community(mst, steps=i)
  #     class <- comm$membership
  #     colors_class <- tableau_color_pal("tableau20")(length(unique(class)))[as.numeric(class)]
  #     set.seed(3)
  #     plot(mst,vertex.size=5,vertex.color=colors_class, vertex.label=NA, main=length(unique(class)))
  #   }
  mat_clus <- matrix(nrow=1000,ncol=1)
  for(i in seq(from=1, to=1000,100)){
    comm <- walktrap.community(mst, steps=i)
    class <- comm$membership
    mat_clus[i,1] <- length(unique(class))
  }
  plot(mat_clus, xlim=c(0,1000))
  #####################################################

  library(gridExtra)
  comm <- walktrap.community(mst, steps=1)
  class <- comm$membership
  names(class) <- comm$names
  colors_class <- tableau_color_pal("tableau20")(length(unique(class)))[as.numeric(class)]
  par(mfcol=c(1,3))
  set.seed(123)
  plot(mst,vertex.size=5,vertex.color=colors_class, vertex.label=NA)
  df <- as.data.frame(table(comm$membership))
  colnames(df) <- c("Cluster", "Cells")
  grid.table(df)
  plot(1,1,axes=F,type="n", xlab="", ylab="")
  set.seed(123)
  plot(mst,vertex.size=5,vertex.color=NA, vertex.label=class)
}

plot.comm <- function(mst,steps){
  library(gridExtra)
  comm <- walktrap.community(mst, steps=steps)
  class <- comm$membership
  names(class) <- comm$names
  colors_class <- tableau_color_pal("tableau20")(length(unique(class)))[as.numeric(class)]
  par(mfcol=c(1,3))
  set.seed(123)
  plot(mst,vertex.size=5,vertex.color=colors_class, vertex.label=NA)
  df <- as.data.frame(table(comm$membership))
  colnames(df) <- c("Cluster", "Cells")
  grid.table(df,col=colors_class)
  plot(1,1,axes=F,type="n", xlab="", ylab="")
  set.seed(123)
  plot(mst,vertex.size=5,vertex.color=NA, vertex.label=class)
}

# data must have cell.name and a column named "order" which
# defines the rank of the cells for the Kolmogorov-Smirnov test
do.gene.ks <- function(data, expr.cutoff=c(2,0)){
  data.ks = copy(data)
  setkey(data.ks,gene,order)
  data.ks[,numcells.on:=sum(log10.cpm>expr.cutoff[1]),by=gene]
  data.ks[,numcells.off:=sum(log10.cpm<=expr.cutoff[2]),by=gene]
  data.ks[,min.numcells:=pmin(numcells.on,numcells.off),by=gene]
  data.ks[,max.numcells:=pmax(numcells.on,numcells.off),by=gene]
  data.ks.2 <- data.ks[min.numcells>0.1*numcells,]
  print(paste("Number of Genes to Analyze: ",length(unique(data.ks.2[,gene]))))

  # Calculate step sizes for each ks site
  data.ks.2[,S_b:=2/numcells]
  data.ks.2[,S_a:=S_b*min.numcells/max.numcells]

  # Calculate steps for each cell and ks site
  # step is 0 if the cell expresses neither or both
  data.ks.2[max.numcells==numcells.on, s_a:=S_a*(log10.cpm>expr.cutoff[1])]
  data.ks.2[max.numcells==numcells.off, s_a:=S_a*(log10.cpm<=expr.cutoff[2])]
  data.ks.2[max.numcells==numcells.on, s_b:=-S_b*(log10.cpm<=expr.cutoff[2])]
  data.ks.2[max.numcells==numcells.off, s_b:=-S_b*(log10.cpm>expr.cutoff[1])]

  # Calculate cumulative distribution function for each ks site
  # Note - it's necessary that data.ks.2 and data.s are in the order specified by "order"!
  setkey(data.ks.2, order, gene)
  data.ks.2[,ks_dist:=as.double(NA)]
  order.g <- unique(data.ks.2[,order])
  print(order.g)

  # Necessary that each cell have a value (either 1 or 0 ) for each gene, or loop will fail
  # and some ks.dist values will be NA
  for(i in 1:length(order.g)){
    if(i==1){
      data.ks.2[order==order.g[i], ks_dist:=s_a+s_b]
    }
    else {
      data.ks.2[order==order.g[i], ks_dist:=data.ks.2[order==order.g[i-1],ks_dist]+s_a+s_b]
    }
  }

  setkey(data.ks.2,gene,order)
  data.ks.2[,d_n_max:=max(ks_dist,na.rm = T),by=gene]
  data.ks.2[,d_n_min:=min(ks_dist,na.rm=T),by=gene]

  data.ks.2
}

# takes output from above and plots d_n_max vs d_n_min with a line
# to separate out the most significant genes
get.best.gene.ks <- function(data.ks, slope,intercept){
  gene.stats <- unique(data.ks[,.(gene,min.numcells,max.numcells,d_n_max,d_n_min,S_b,S_a)])

  plt1=ggplot(gene.stats,aes(x=d_n_max,y=d_n_min,color=log2(min.numcells)))+geom_point()+geom_abline(intercept=intercept, slope=slope)+
    theme_bw()
  print(plt1)
  gene.stats[,is_best:=(intercept+slope*d_n_max)>d_n_min]
  plt2=ggplot(gene.stats,aes(x=d_n_max,y=d_n_min,color=is_best))+geom_point()+geom_abline(intercept=intercept, slope=slope)+
    theme_bw()

  data.best <- data.ks[gene %in% gene.stats[is_best==T,gene]]
  list(data.best, plt1,plt2)
}


rescale <- function(x){
  x <- (x-min(x))/(max(x)-min(x))
  x
}

shorten_cellname <- function(cell.name){
  cell.name_short:=gsub("_","-",cell.name)
  cell.name_short:=gsub("-IL.*","",cell.name_short)
  cell.name_short
}

legend.col <- function(col, lev){

  opar <- par

  n <- length(col)

  bx <- par("usr")

  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n

  xx <- rep(box.cx, each = 2)

  par(xpd = TRUE)
  for(i in 1:n){

    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])

  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}

do.deseq <- function(cast.group1, cast.group2, filename){
  cast.group1 <- data.table(cast.group1,keep.rownames = T)
  colnames(cast.group1)[colnames(cast.group1)=="rn"] <- "cell.name"
  cast.group2 <- data.table(cast.group2,keep.rownames = T)
  colnames(cast.group2)[colnames(cast.group2)=="rn"] <- "cell.name"
  cast.group1[, group:="group1"]; cast.group2[, group:="group2"]

  cast.all <- rbindlist(list(cast.group1, cast.group2),use.names = T,fill=F)

  cast.all <- cast.all[, c("cell.name","group",colnames(cast.all)[!colnames(cast.all) %in% c("cell.name","group")]),with=F] # move group to front

  counts.all <- as.data.frame(cast.all); rownames(counts.all)<-counts.all$cell.name
  counts.all <- counts.all[,3:ncol(counts.all)]
  counts.all <- as.data.frame(t(counts.all))
  tmp <- sapply(counts.all, FUN = as.integer); rownames(tmp) <- rownames(counts.all)
  counts.all <- tmp; rm(tmp)

  library(DESeq2)
  colData=as.data.frame(cast.all[, group], row.names=cast.all[,cell.name]);
  colnames(colData)<-"group"
  dds=DESeq2::DESeqDataSetFromMatrix(countData=counts.all,
                                     colData=colData,
                                     design= ~group)

  gc()
  dds <- DESeq2::DESeq(dds)
  res.all <- DESeq2::results(dds)
  res.all <- res.all[order(res.all$padj),]
  res.all <- res.all[!is.na(res.all$log2FoldChange), ]

  write.csv(res.all, file=filename, quote=F)

  res.all
}

# function to return just the p-value of the Kruskal-Wallis test
kw.pval = function(x, y){
  # res = kruskal.test(list(x, y))
  res = t.test(x, y)
  res$p.value
}

# do differential expression comparing 2 groups of cells with Kruskal-Wallis test
# Returns a list - first element is the differential expression, second element is a plot
diff.expr.kwtest <- function(data, cells.1, cells.2, n.cores = 8, from.data.table = F){
  library(parallel)
  library(pbapply)
  library(ggrepel)

  if(from.data.table){
    mat = cast.to.df(data[cell.name %in% c(cells.1, cells.2)], to.matrix = T)
  } else{
    mat <- data[unique(c(cells.1, cells.2)), ]
  }

  # For non-detected genes, need to gracefully return an NA value
  genes.all = colnames(mat)
  genes.use = colnames(mat)[colSums(mat > 0) > 0]
  genes.exclude = genes.all[!genes.all %in% genes.use]

  genes.use = colnames(mat)
  itr.fxn = mclapply

  p.val = unlist(itr.fxn(genes.use, function(x) kw.pval(mat[cells.1, x], mat[cells.2, x]), mc.cores = n.cores))

  to.return=data.table(p.val,gene = genes.use)
  to.return[, p.adj := p.adjust(p.val, method = 'fdr')]
  setorder(to.return, p.adj)


  diff=unlist(itr.fxn(genes.use, function(x) mean(mat[cells.1, x]) - mean(mat[cells.2, x]), mc.cores = n.cores))
  diff.return=data.table(diff, gene = genes.use)

  len.1 = length(cells.1)
  len.2 = length(cells.2)
  diff.frxn=unlist(itr.fxn(genes.use, function(x) sum(mat[cells.1, x] > 0)/len.1 - sum(mat[cells.2, x] > 0)/len.2,
                           mc.cores = n.cores))
  diff.frxn.return=data.table(diff.frxn, gene = genes.use)

  to.return = merge(to.return, diff.return, by = 'gene')
  to.return = merge(to.return, diff.frxn.return,   by = 'gene')

  setorder(to.return, diff)
  plt = ggplot(to.return, aes(diff, log10(p.adj))) + geom_point() +
    theme_bw() +
    # geom_point(aes(diff, log10(p.adj)), data = to.return[abs(diff) > .3][p.adj < 0.05]) +
    xlab('Expr. Diff. (log10)') +
    # ylab('p (KW-test, BY adjusted)') +
    ylab('log10 fdr') +
    geom_label_repel(aes(diff, log10(p.adj), label = gene), data = to.return[order(p.adj)][1:10])


  # Add the non-expressed genes back into the data.
  # Give values of 0 for differential expr and NA for pvalues
  # Make sure result dataframe has same gene order as input matrix
  if(length(genes.exclude) > 0){
    tmp = data.frame(gene = genes.exclude, p.val = NA, p.adj = NA, diff = 0, diff.frxn = 0)
    to.return = rbind(to.return, tmp[, colnames(to.return)])
    to.return = to.return[match(genes.all, to.return$gene), ]

  }

  return(list(to.return, plt))

}

# function to find unique genes in one population relative to N others
# UNFINISHED - work in progress
find.unique.markers = function(data, cells.target, cells.rest.list, n.cores = 8){

  g.use = unique(as.character(data[, gene]))
  # differential comparison with each population
  i = 0
  for(cells.rest.i in cells.rest.list){
    res.i = diff.expr.kwtest(data, cells.target, cells.rest.i, n.cores = n.cores)
    res.i[[1]] = res.i[[1]][!is.na(p.adj)]
    g.use = g.use[g.use %in% res.i[[1]][, gene]]
    if(i == 0){
      diff.list = list(res.i[[1]])
    } else {
      diff.list = list(diff.list, res.i[[1]])
    }
    i = i+1
  }

  # get minimum significance and minimum differential expression
  # first, upregulated genes
  # only look at genes w/non-NA values
  for(g in g.use){

  }


}


cor.pval <- function(x,y,method='spearman'){
  cor.test(x,y,method=method)$p.value
}

plot.ggpairs.wExptCols <- function(cell.info, dir='',num.pcs=6,suffix='',height=21,width=21,suffix.2='',sigType='PC',annot.expt=F,size=1,
                                   plt.pos.neg=F, plot.read.depth=F, feature.annot = 'experiment'){


  # make sure directory name ends in slash (need to change to run on Windows!)
  if(substr(dir,nchar(dir),nchar(dir))!="/"){
    dir <- paste0(dir,"/")
  }
  library(GGally)
  library(car)
  setkey(data)

  if(sigType=='PC'){
    pc.dims <- paste0('PC',rep(1:num.pcs))
    pc.scores <- fread(paste0(dir,'PC_allscores',suffix,'.csv'));
    pc.scores <- pc.scores[,c('cell.name',pc.dims),with=F]
  }
  else if(sigType=='Sig'){
    if(plt.pos.neg){
      pc.dims <- c(paste0('PC',rep(1:num.pcs,times=rep(2,num.pcs)),c('.pos','.neg')))
      pc.scores <- fread(paste0(dir,'PC_allscores',suffix,'.csv'));
      pc.scores <- pc.scores[,c('cell.name',pc.dims),with=F]

    }
    else{
      pc.dims <- paste0('PC',rep(1:num.pcs),'.score')
      pc.scores <- fread(paste0(dir,'PC_allscores',suffix,'.csv'));
      pc.scores <- pc.scores[,c('cell.name',pc.dims),with=F]; setkey(pc.scores)

    }
  }
    setkey(pc.scores)


  if(plot.read.depth){
    pc.scores = merge(pc.scores, cell.info[, .(cell.name, sum.counts)], by = 'cell.name')
    pc.scores[, read.depth := log10(sum.counts + 1)]
    pc.scores[, sum.counts:=NULL]
  }

  if(annot.expt) {
    cell.info[, experiment := get(feature.annot)]
    pc.scores <- merge(pc.scores,cell.info[,.(cell.name,num.genes,experiment)],by='cell.name')
    setorder(pc.scores,cell.name)
    plt.data = as.data.frame(pc.scores[, -c('cell.name','experiment'), with = F])
    print(colnames(plt.data))
    plt.annot <- as.character(pc.scores[['experiment']])
    plt.color <- gg_color_hue(length(unique(plt.annot)))[as.factor(plt.annot)]

    jpeg(file=paste0(dir,"ggpair_",sigType,suffix,suffix.2,".jpg"),height=height,width=width,units='in',res=250)
    # print(ggpairs(pc.scores[,c(pc.dims[1:num.pcs],'num.genes','experiment'),with=F],mapping = aes(color=experiment))+theme_bw())
    print(plot(plt.data,pch=20,cex=size,col=plt.color))#,xlim=c(0,1),ylim=c(0,1)))
    # print(spm(pc.scores[,c(pc.dims,'num.genes'),with=F],pch=20,cex=size,groups=pc.scores[,experiment],
    #           smooth=F,reg.line=F))
    dev.off()
  }
  else {
    pc.scores <- merge(pc.scores,cell.info[,.(cell.name,num.genes)],by='cell.name')
    setorder(pc.scores,cell.name)
    plt.data = as.data.frame(pc.scores[, colnames(pc.scores) != 'cell.name', with = F])
    print(colnames(plt.data))

    jpeg(file=paste0(dir,"ggpair_",sigType,suffix,suffix.2,".jpg"),height=height,width=width,units='in',res=250)
    # print(ggpairs(pc.scores[,c(pc.dims[1:num.pcs],'num.genes'),with=F])+theme_bw())
    print(plot(plt.data,pch=20,cex=size,col='black'))#,xlim=c(0,1),ylim=c(0,1)))
    dev.off()
  }



}

# make biplots of gene expression signatures (librarys two character vectors of genes)
sig.biplot <- function(data,fname,genes.pos,genes.neg){
  # generate signature value of top PC +/- correlated genes
  genes.use <- c(genes.pos,genes.neg)
  data.sig <- data[gene %in% genes.use];
  data.sig[gene %in% genes.pos,genetype:="pos"];data.sig[gene %in% genes.neg,genetype:="neg"]
  data.sig[,numgenes:=length(unique(gene)),by=.(genetype)]
  data.sig[,scaled_exp:=log10.cpm/max(log10.cpm)/numgenes,by=.(gene)]
  data.sig[,pc.sig:=sum(scaled_exp),by=.(cell.name,genetype)];

  cast.sig <- dcast.data.table(data.sig, cell.name +experiment + experiment_color ~ genetype, fun.aggregate = mean, value.var = "pc.sig",fill = 0)

  # make color scale
  setkey(data.sig)
  expt.cols <- unique(data.sig[,.(experiment,experiment_color)]);
  expt2col <- expt.cols[['experiment_color']]; names(expt2col)<-expt.cols[['experiment']]

  plt=ggplot(cast.sig,aes(x=pos,y=neg,color=experiment))+geom_point()+theme_bw()+
    scale_color_manual(values=expt2col)
  ggsave(plt,file=fname,height=4,width=6)
}

sig.histogram <- function(data, fname,genes.pos,genes.neg,num.bins){
  # generate signature value of top PC +/- correlated genes
  cast.sig <- get.sig(data,genes.pos=genes.pos,genes.neg=genes.neg,write_csv = F)
  cast.sig[,frac.pos:=pos/(pos+neg)]
  # plot histograms of signature values
  range = max(cast.sig[['frac.pos']])-min(cast.sig[['frac.pos']])
  plt=ggplot(cast.sig,aes(x=frac.pos,y=..count..))+geom_histogram(binwidth=range/num.bins)+theme_classic()
  ggsave(plt,file=fname,height=4,width=5)

}

get.sig <- function(data,genes.pos,genes.neg,write_csv=F,fname='',scale=T,center=F,normalize=F, from.mat = T,
                    dt.out = F){
  # generate signature value of top PC +/- correlated genes
  genes.use <- c(genes.pos,genes.neg)

  if(!from.mat){
    setkey(data)
    # get all cell names to avoid cells being removed that don't express any of genes.pos or genes.neg
    cellnames=unique(data[,.(cell.name)])
    cells.noExpr=cellnames[!cell.name %in% unique(data[gene %in% genes.use,cell.name]),cell.name]
    data[cell.name %in% cells.noExpr][gene==genes.pos[1],log10.cpm:=0]
    data[cell.name %in% cells.noExpr][gene==genes.neg[1],log10.cpm:=0]

    data.sig <- data[gene %in% genes.use,.(cell.name,gene,log10.cpm)];setkey(data.sig)
    cast.melt <- dcast.data.table(data.sig,cell.name~gene,value.var = 'log10.cpm',fill=0)
    data.sig <- melt(cast.melt,value.name = 'log10.cpm',id.vars = 'cell.name',variable.name = 'gene')
    setkey(data.sig)
    data.sig[gene %in% genes.pos,genetype:="pos"];data.sig[gene %in% genes.neg,genetype:="neg"]

    data.sig[,numgenes:=length(unique(gene)),by=.(genetype)]
    data.sig[,scaled_exp:=log10.cpm/max(log10.cpm),by=.(gene)]
    data.sig[,centered_exp:=log10.cpm-mean(log10.cpm),by=.(gene)]
    data.sig[,centered_scaled_exp:=scaled_exp-mean(scaled_exp),by=.(gene)]
    data.sig[,norm_exp:=scale(log10.cpm,center = T,scale = T),by=.(gene)]

    if(normalize){data.sig[,pc.sig:=sum(norm_exp),by=.(cell.name,genetype)] }
    else if(scale && !center) { data.sig[,pc.sig:=sum(scaled_exp)/numgenes,by=.(cell.name,genetype)] }
    else if(center && !scale){data.sig[,pc.sig:=sum(centered_exp)/numgenes,by=.(cell.name,genetype)]}
    else if(center && scale){data.sig[,pc.sig:=sum(centered_scaled_exp),by=.(cell.name,genetype)]}
    else{data.sig[,pc.sig:=sum(centered_exp),by=.(cell.name,genetype)]}

    cast.sig <- dcast.data.table(data.sig, cell.name ~ genetype, fun.aggregate = mean, value.var = "pc.sig",fill = 0)
    setkey(cast.sig)

  }
  else{
    data.scale = data[, genes.use] / apply(data[, genes.use], 2, max)

    if(dt.out){
      cast.sig = data.table(cell.name = rownames(data.scale),
                            neg = rowSums(data.scale[, genes.neg]),
                            pos = rowSums(data.scale[, genes.pos]))
      cast.sig[, pos := pos/max(pos)]
      cast.sig[, neg := neg/max(neg)]
      setkey(cast.sig)
    }
    else { # as of 170723 default behavior is to output a data.frame
      cast.sig <- data.frame(neg = rowSums(data.scale[, genes.neg]),
                            pos = rowSums(data.scale[, genes.pos]))
      rownames(cast.sig) <- rownames(data.scale)
      cast.sig$pos <- cast.sig$pos / max(cast.sig$pos)
      cast.sig$neg <- cast.sig$neg / max(cast.sig$neg)
    }

  }


  if(write_csv){
    write.csv(cast.sig,fname,quote=F,row.names=T)
  }

  cast.sig

}

get.2dsig <- function(data,pcs.use=c('PC1','PC2'),dir.pc,num.genes=30,write_csv=F,fname=''){
  pc.loadings <- fread(paste0(dir.pc,'PCloadings.csv'))

  setorderv(pc.loadings,pcs.use[1]); genes.pc1.pos <- pc.loadings[,gene][1:num.genes]
  setorderv(pc.loadings,pcs.use[1],order = -1); genes.pc1.neg <- pc.loadings[,gene][1:num.genes]
  setorderv(pc.loadings,pcs.use[2]); genes.pc2.pos <- pc.loadings[,gene][1:num.genes]
  setorderv(pc.loadings,pcs.use[2],order = -1); genes.pc2.neg <- pc.loadings[,gene][1:num.genes]

  sig.pc1 <- get.sig(data,genes.pc1.pos,genes.pc1.neg)
  sig.pc2 <- get.sig(data,genes.pc2.pos,genes.pc2.neg)
  sig.pc1[,pc1.score:=pos-neg]; sig.pc2[,pc2.score:=pos-neg]
  sig.scores <- merge(sig.pc1[,.(cell.name,pc1.score)],sig.pc2[,.(cell.name,pc2.score)],by='cell.name')
  setnames(sig.scores,c('cell.name',paste0(pcs.use[1],'.score'),paste0(pcs.use[2],'.score')))

  sig.scores

}

# calculate a signature for any set of cells from the genes correlated to any PC
get.sig.fromPC <- function(data, pc.dir, pcs.pos, pcs.neg, genes.per.pc = 30, remove.shared.genes = T){
  library(data.table)
  library(plyr)
  loadings = fread(file.path(pc.dir, 'PCloadings.csv'))

  # still need to make this compatible with PC.score = PC.pos - PC.neg
  genes.pos <- get.pc.genes(loadings, pcs.pos, ng = genes.per.pc)
  genes.neg <- get.pc.genes(loadings, pcs.neg, ng = genes.per.pc)

  if(remove.shared.genes){
    genes.shared = intersect(genes.pos, genes.neg)
    genes.pos  = genes.pos[!genes.pos %in% genes.shared]
    genes.neg  = genes.neg[!genes.neg %in% genes.shared]
  }

  # calculate signature
  genes.use <- unique(c(genes.pos,genes.neg))
  data.g.use = data[gene %in% genes.use]
  cast.sig <- get.sig(data.g.use, genes.pos, genes.neg)
  cast.sig[, score := pos - neg]
  cast.sig = cast.sig[, .(cell.name, pos, neg, score)]
  return(cast.sig)
}


# Usage
#gsigsave(cast.plot,'PC1.pos','PC1.neg','Drd1a')
gsigsave <- function(data,x,y,gene,size=3.5,palette="YlOrRd",height=4,width=6,dir.use=NA,plt.legend=F,plt.axis.titles=T) {
  rm(plt)


  if(!is.na(dir.use)){
    cordir=dir.use
  }
  d.plot <- data.table(x.plot=data[[x]],y.plot=data[[y]],g.plot=data[[gene]])
  plt=ggplot(d.plot,aes(x.plot,y.plot,color=g.plot))+ geom_point(size=size,pch=20)+
    scale_color_gradientn(colours=brewer.pal(9,palette)[4:9])+xlab(x)+ylab(y)+theme_classic()+
    theme(axis.title=element_text(size=18),axis.text=element_text(size=18))
  if(!plt.legend){
    plt=plt+theme(legend.position='none')
  }
  if(!plt.axis.titles){
    plt=plt+theme(axis.title=element_blank())
  }

  ggsave(plt,file=paste0(cordir,x,'_',y,'-',gene,'.pdf'),height=height,width=width)
}

gsigsave_spec <- function(data,x,y,gene,size=3.5,palette="Spectral",height=4,width=6,dir.use=NA) {
  rm(plt)

  if(!is.na(dir.use)){
    cordir=dir.use
  }
  d.plot <- data.table(x.plot=data[[x]],y.plot=data[[y]],g.plot=data[[gene]])
  plt=ggplot(d.plot,aes(x.plot,y.plot,color=g.plot))+ geom_point(size=size,pch=20)+
    scale_color_gradientn(colours=brewer.pal(11,palette)[11:2])+xlab(x)+ylab(y)+theme_classic()+
    theme(axis.title=element_text(size=18),axis.text=element_text(size=14),legend.title=element_blank())

  ggsave(plt,file=paste0(cordir,x,'_',y,'-',gene,'.pdf'),height=height,width=width)
}

make.sigbiplots <- function(data,pc.use,pc.dir,size=3.5,height=4,width=6,num.genes=30,pc.loadings.fname='PCloadings.csv') {
  dir.use=paste0(pc.dir,'biplots/')

  pc.loadings <- fread(paste0(pc.dir,pc.loadings.fname))
  setorderv(pc.loadings,pc.use); genes.pos <- pc.loadings[,gene][1:num.genes]
  setorderv(pc.loadings,pc.use,order=-1); genes.neg <- pc.loadings[,gene][1:num.genes]
  sig.pc <- get.sig(data,genes.pos,genes.neg)
  setnames(sig.pc,c('cell.name',paste0(pc.use,'.neg'),paste0(pc.use,'.pos')))

  genes.plot=c(genes.pos,genes.neg,'num.genes')
  cast.plot <- dcast.data.table(data[gene %in% genes.plot],cell.name+num.genes+experiment~gene,value.var='log10.cpm',fill = 0);cast.plot <- merge(cast.plot,sig.pc,by='cell.name')


  for(g in genes.plot){
    gsigsave(cast.plot,paste0(pc.use,'.pos'),paste0(pc.use,'.neg'),g,size=size,height=height,width=width,dir.use = dir.use)
  }
  cast.plot <- dcast.data.table(data,cell.name+num.genes+experiment~gene,value.var='log10.cpm',fill = 0);cast.plot <- merge(cast.plot,sig.pc,by='cell.name')
  cast.plot
}

# data is a matrix of log10 expression OR a matrix/dataframe with a categorical/discrete annotation
# and row.names = cell name
pc.plot <- function(data, annot, pc.use, pc.dir, suffix = '', is.gene = T){
  pc.scores <- read.csv(paste0(pc.dir,'PC_allscores',suffix,'.csv'))
  data = data[rownames(data) %in% as.character(pc.scores[, 'cell.name']), ]
  g.expr = merge(data.frame(cell.name = rownames(data), annot = data[, annot]),
                 pc.scores[, c('cell.name', pc.use)],
                 by = 'cell.name')
  colnames(g.expr) <-  c('cell.name', 'annot', pc.use)
  # 3D plot
  if(length(pc.use) == 3){
    library(plot3Drgl)

    if(is.gene){
      plot3d(g.expr[, pc.use[1]], g.expr[, pc.use[2]], g.expr[, pc.use[3]])#, col = vecs2rgb(annot)))

    } else {
      g.expr[, annot] <- gg_color_hue(length(unique(g.expr[, annot])))[factor(g.expr[, annot])]
      with(g.expr, plot3d(`pc.use[1]`, `pc.use[2]`, `pc.use[3]`, col = `annot`))

    }

  } else{ #2D plot
    if(is.gene){
      p = ggplot(g.expr, aes_string(pc.use[1], pc.use[2])) +
        geom_point(size=1, shape = 1) +
        geom_point(aes(color=annot), size=1, shape = 16) +
        scale_color_gradientn(colours=c(brewer.pal(9,'YlOrBr')[c(2:9)]), guide = F)

    } else {
      p = ggplot(g.expr, aes_string(pc.use[1], pc.use[2], color = annot)) +
        geom_point(size=1, shape = 1) +
        geom_point(aes(color=annot), size=1, shape = 16) +
        theme(legend.position = 'none')

    }

    p
  }

}


make.facet.biplots <- function(data,pc.use,pc.dir,size=3.5,height=12,width=16,num.genes=8,pc.loadings.fname='PCloadings',pc.scores.fname='PC_allscores',
                               method.clust='pearson',ax.suf=c('.score','.score'),scores.dir=NA,print.dir=NA,scores.pc='',scale=F,custom.genes='',
                               suffix='', plot.pdf = T, fixedcoord = T, from.mat = F) {

  if(length(pc.use) == 1) {
    pc.use = c(pc.use, pc.use)
    ax.suf <- c('.pos','.neg')
  }
  if(is.na(scores.dir)){
    scores.dir=pc.dir
    scores.pc=pc.use
  }
  if(is.na(print.dir)){
    print.dir=paste0(pc.dir,'biplots/')
  }

  if(custom.genes[1]==''){
    pc.loadings <- fread(paste0(pc.dir,pc.loadings.fname,suffix,'.csv'))

    setorderv(pc.loadings,pc.use[1]); genes.pc1.pos <- pc.loadings[,gene][1:num.genes]
    setorderv(pc.loadings,pc.use[1],order = -1); genes.pc1.neg <- pc.loadings[,gene][1:num.genes]
    genes.plot <- c(genes.pc1.pos,genes.pc1.neg)

    if(pc.use[1]!=pc.use[2]){
      setorderv(pc.loadings,pc.use[2]); genes.pc2.pos <- pc.loadings[,gene][1:num.genes]
      setorderv(pc.loadings,pc.use[2],order = -1); genes.pc2.neg <- pc.loadings[,gene][1:num.genes]
      genes.plot <- c(genes.pc1.pos,genes.pc1.neg,genes.pc2.pos,genes.pc2.neg)
    }
  }

  pc.scores <- fread(paste0(scores.dir,pc.scores.fname,suffix,'.csv'))

  if(custom.genes[1]=='') {
    if(from.mat) data <- as.data.table(melt(data[, unique(genes.plot)], varnames = c('cell.name', 'gene'),
                              value.name = 'log10.cpm'))
    if(class(data)[1] != 'data.table') data <- as.data.table(data)
    data.plot <- merge(data[gene %in% genes.plot], pc.scores, by='cell.name')


  }
  else {
    custom.genes=as.character(custom.genes)
    if(from.mat) data <- as.data.table(melt(data[, unique(custom.genes)], varnames = c('cell.name', 'gene'),
                                         value.name = 'log10.cpm'))
    if(class(data)[1] != 'data.table') data <- as.data.table(data)
    data[, gene:=as.character(gene)]
    data.plot <- merge(data[gene %in% custom.genes],pc.scores,by='cell.name')
  }

  # only plot genes expressed in > 0 cells! Otherwise spearman correlation throws an error
  data.plot[, numcells := sum(log10.cpm > 0), by = 'gene']
  data.plot <- data.plot[numcells > 0]

  if(length(unique(data.plot[,gene]))>2) {
    data.plot[, gene:=factor(gene, levels = get.clustorder(data.plot,method.cor=method.clust))]
  }

  if(scale){
    data.plot[,log10.cpm:=(log10.cpm-min(log10.cpm))/(max(log10.cpm)-min(log10.cpm)),by=gene]
  }
  plt=ggplot(data.plot,aes_string(paste0(scores.pc[1],ax.suf[1]),paste0(scores.pc[2],ax.suf[2])))+
    # geom_point(size=size, shape = 1) +
    geom_point(aes(color=log10.cpm), size=size, shape = 16) +
    theme_bw() +
    facet_wrap(~gene)+
    theme(strip.text=element_text(size=20),axis.title=element_text(size=20),legend.text=element_text(size=18),
          panel.grid = element_blank()) +
    # scale_color_gradientn(colours=c(brewer.pal(9,'YlOrRd')[c(2:8,8,9,9)],'black','black'))
    # scale_color_gradientn(colours=c(brewer.pal(9,'YlOrBr')[c(2:9)]))
    scale_color_gradientn(colours=brewer.pal(9,'YlOrBr')[c(3:9)])
    # scale_color_gradientn(colours=c(brewer.pal(9,'YlOrBr')[c(3:9)],'black'))

  if(fixedcoord) plt = plt + coord_fixed()

  # if(length(custom.genes) > 0) {ggsave(plt,filename = paste0(print.dir,'facet_biplot_',pc.use[1],ax.suf[1],'_',pc.use[2],ax.suf[2],'_',custom.genes[1],'.pdf'),height=height,width=width)}
  # else{ggsave(plt,filename = paste0(print.dir,'facet_biplot_',pc.use[1],ax.suf[1],'_',pc.use[2],ax.suf[2],'.pdf'),height=height,width=width)}
  if(plot.pdf == T){
    if(custom.genes[1]=='') {ggsave(plt,filename = paste0(print.dir,'facet_biplot_',pc.use[1],ax.suf[1],'_',pc.use[2],
                                                          ax.suf[2],suffix,'_',custom.genes[1],'.pdf'),height=height,width=width)} else{
                                                            ggsave(plt,filename = paste0(print.dir,'facet_biplot_',pc.use[1],ax.suf[1],'_',pc.use[2],
                                      ax.suf[2],custom.genes[1],suffix,'.pdf'),height=height,width=width)}
  } else{
    if(custom.genes[1]=='') {ggsave(plt,filename = paste0(print.dir,'facet_biplot_',pc.use[1],ax.suf[1],'_',pc.use[2],ax.suf[2],suffix,'_',custom.genes[1],'.jpg'),
                                    height=height,width=width,dpi = 250)} else{
                                      ggsave(plt,filename = paste0(print.dir,'facet_biplot_',pc.use[1],ax.suf[1],'_',
                                                                   pc.use[2],ax.suf[2],custom.genes[1],suffix,'.jpg'),height=height,width=width,dpi=250)}

  }


}

# Calculate two dimensions based on top N genes from PC loadings from different analyses (specified by directory)
# and print to directory of x dimension unless specified.
biplots.mixed.scores <- function(data,pc.use,pc.dir,size=3.5,height=28,width=34,num.genes=40,pc.loadings.fname='PCloadings.csv',pc.scores.fname='PC_allscores.csv',
                               method.clust='pearson',ax.suf=c('.score','.score'),print.dir=NA,scale=F,custom.genes='') {
  scores.dir=pc.dir
  scores.pc=pc.use

  if(is.na(print.dir)){
    print.dir=paste0(pc.dir[1],'biplots/')
  }

  pc.loadings.1 <- fread(paste0(pc.dir[1],pc.loadings.fname))
  setorderv(pc.loadings.1,pc.use[1]); genes.pc1.pos <- pc.loadings.1[,gene][1:num.genes]
  setorderv(pc.loadings.1,pc.use[1],order = -1); genes.pc1.neg <- pc.loadings.1[,gene][1:num.genes]
  sig.1=get.sig(data,genes.pc1.pos,genes.pc1.neg)
  if(ax.suf[1]=='.score') {sig.1[,sig.1.score:=pos-neg]}
  else if(ax.suf[1]=='.pos') {sig.1[,sig.1.score:=pos]}
  else if(ax.suf[1]=='.neg') {sig.1[,sig.1.score:=neg]}

  pc.loadings.2 <- fread(paste0(pc.dir[2],pc.loadings.fname))
  setorderv(pc.loadings.2,pc.use[2]); genes.pc2.pos <- pc.loadings.2[,gene][1:num.genes]
  setorderv(pc.loadings.2,pc.use[2],order = -1); genes.pc2.neg <- pc.loadings.2[,gene][1:num.genes]
  sig.2=get.sig(data,genes.pc2.pos,genes.pc2.neg)
  if(ax.suf[2]=='.score') {sig.2[,sig.2.score:=pos-neg]}
  else if(ax.suf[2]=='.pos') {sig.2[,sig.2.score:=pos]}
  else if(ax.suf[2]=='.neg') {sig.2[,sig.2.score:=neg]}

  pc.scores <- merge(sig.1[,.(cell.name,sig.1.score)],sig.2[,.(cell.name,sig.2.score)],by='cell.name')
  genes.plot <- c(genes.pc1.pos,genes.pc1.neg,genes.pc2.pos,genes.pc2.neg)

  if(custom.genes[1]=='') {data.plot <- merge(data[gene %in% genes.plot],pc.scores,by='cell.name')}
  else {data.plot <- merge(data[gene %in% custom.genes],pc.scores,by='cell.name')}

  if(length(unique(data.plot[,gene]))>2) {
    data.plot[,gene:=factor(gene,levels=get.clustorder(data.plot,method.cor=method.clust))]
  }

  if(scale){ data.plot[,log10.cpm:=(log10.cpm-min(log10.cpm))/(max(log10.cpm)-min(log10.cpm)),by=gene] }

  plt=ggplot(data.plot,aes(sig.1.score,sig.2.score,color=log10.cpm))+geom_point(size=size)+theme_classic()+
    scale_color_gradientn(colours=brewer.pal(9,'YlOrBr')[4:9])+facet_wrap(~gene)+theme(strip.text=element_text(size=20),axis.title=element_text(size=20))

  # remove trailing "/" in directory names
  for(i in 1:length(pc.dir)) { pc.dir[i]=substring(pc.dir[i],1,nchar(pc.dir[i])-1) }

  if(length(custom.genes) > 0) {
    ggsave(plt,filename = paste0(print.dir,'mixed_',pc.dir[1],'-',pc.use[1],ax.suf[1],'_',
                                 pc.dir[2],'-',pc.use[2],ax.suf[2],'_',custom.genes[1],'.pdf'),height=height,width=width)
    }
  else{
    ggsave(plt,filename = paste0(print.dir,'mixed_',pc.dir[1],'-',pc.use[1],ax.suf[1],'_',
                                    pc.dir[2],'-',pc.use[2],ax.suf[2],'.pdf'),height=height,width=width)
    }
}



plot.rmeans <- function(data,sig.score,x,genes.use,k=20,fname='',height=10,width=12,size=1,line=T) {
  cast.plot <- dcast.data.table(data[gene %in% genes.use],cell.name~gene,value.var='log10.cpm',fill = 0);cast.plot <- merge(cast.plot,sig.score,by='cell.name')
  setorderv(cast.plot,x)
  rmeans=as.data.table(rollmean(cast.plot[,2:ncol(cast.plot),with=F],k = k))
  # rsd=as.data.table(rollapply(cast.plot[,2:ncol(cast.plot),with=F],k =k))
  d.rmeans <- melt(rmeans,id.vars = x,variable.name = 'gene',value.name = 'log10.cpm')
  d.rmeans[,gene:=factor(gene,levels=genes.use)]
  plt=ggplot(d.rmeans,aes_string(x,'log10.cpm'))+geom_point(size=size)+theme_classic()+facet_wrap(~gene)
  if(line){
    plt=plt=geom_line()
  }
  ggsave(plt,file=fname,height=height,width=width)
}

rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE,
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{

  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)

  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0),
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)

  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)

  # Add plane
  if(show.plane)
    xlim <- xlim/1.1; zlim <- zlim /1.1
  rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
             z = c(zlim[1], zlim[2], zlim[2], zlim[1]))

  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5,
             emission=bbox.col[1], specular=bbox.col[1], shininess=5,
             xlen = 3, ylen = 3, zlen = 3)
  }
}

# uses melted data.table or cell.name X gene matrix (set from.matrix to T if using matrix as input)
get.clustorder <- function(data,method='ward',method.cor='pearson',id.var='cell.name',value.name='log10.cpm',
                           from.matrix = F) {

  if(!from.matrix){ # assuming this is from a data.table with cell.name, gene, and expression columns

    data.use <- data[,c(id.var,'gene',value.name),with=F]
    setnames(data.use,c(id.var,'gene',value.name))
    cast.data = as.data.frame(dcast.data.table(data.use, cell.name ~ gene, value.var=value.name,fill=0))
    rownames(cast.data)=cast.data[,id.var];
    cast.data =as.matrix(cast.data[,2:ncol(cast.data)])

  } else{ cast.data = copy(data)}

  # only use genes w/ >0 expression! or error in correlation
  cast.data = cast.data[, apply(cast.data, 2, sum) > 0]
  coldist=as.dist(1-cor(cast.data,method=method.cor));
  hc <- hclust(coldist, method=method)

  library(cba)
  optimal.col <- order.optimal(coldist,hc$merge)
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order
  ddc <- as.dendrogram(hc)
  labels(ddc)
}



get.cellorder <- function(data,method='ward',method.cor='pearson',id.var='cell.name',value.name='log10.cpm') {
  data.use <- data[,c(id.var,'gene',value.name),with=F]
  setnames(data.use,c('cell.name','gene',value.name))
  cast.data = as.data.frame(dcast.data.table(data.use, cell.name ~ gene, value.var=value.name,fill=0))
  rownames(cast.data)=cast.data[,'cell.name'];cast.data=cast.data[,2:ncol(cast.data)]
  coldist=as.dist(1-cor(cast.data,method=method.cor));
  hc <- hclust(coldist, method=method)

  library(cba)
  optimal.col <- order.optimal(coldist,hc$merge)
  hc$merge <- optimal.col$merge; hc$order <- optimal.col$order
  ddc <- as.dendrogram(hc)
  labels(ddc)
}

make.expt.biplot <- function(data,pc.use,pc.dir,size=3.5,height=3.5,width=5,pc.scores.fname='PC_allscores',
                              ax.suf=c('.score','.score'),scores.dir=NA,print.dir=NA,scores.pc='',suffix='') {
  if(is.na(scores.dir)){
    scores.dir=pc.dir
    scores.pc=pc.use
  }
  if(is.na(print.dir)){
    print.dir=pc.dir
  }

  cell.info = unique(data[,.(cell.name,experiment)])
  pc.scores <- fread(paste0(scores.dir,pc.scores.fname,suffix,'.csv'))
  cell.info <- merge(cell.info,pc.scores,by='cell.name')
  plt=ggplot(cell.info,aes_string(paste0(scores.pc[1],ax.suf[1]),paste0(scores.pc[2],ax.suf[2]),color='experiment'))+geom_point(size=size)+theme_classic()+
    scale_color_manual(values=expt2col)+theme(axis.text=element_text(size=12))
  ggsave(plt,filename = paste0(print.dir,'biplots/expt_biplot_',pc.use[1],ax.suf[1],'_',pc.use[2],ax.suf[2],'.pdf'),height=height,width=width)
}
make.cond.biplots <- function(data,pc.use,pc.dir,size=3.5,height=6,width=8,pc.scores.fname='PC_allscores.csv',
                              ax.suf=c('.score','.score'),scores.dir=NA,print.dir=NA,scores.pc='') {
  if(is.na(scores.dir)){
    scores.dir=pc.dir
    scores.pc=pc.use
  }
  if(is.na(print.dir)){
    print.dir=pc.dir
  }

  cell.info = unique(data[,.(cell.name,experiment)])
  cell.info = merge(cell.info, expt.cond,by='experiment')
  cell.info[,is_cond:=1]
  cell.info = dcast.data.table(cell.info,cell.name ~ condition, value.var = 'is_cond',fill=0)
  cell.info = melt(cell.info,id.vars=c('cell.name'),variable.name = 'condition',value.name='is_cond')
  pc.scores <- fread(paste0(scores.dir,pc.scores.fname))
  data.plot <- merge(cell.info,pc.scores,by='cell.name')

  plt=ggplot(data.plot,aes_string(paste0(scores.pc[1],ax.suf[1]),paste0(scores.pc[2],ax.suf[2]),color='is_cond'))+geom_point(size=size)+theme_classic()+
    scale_color_gradientn(colours=brewer.pal(9,'YlOrBr')[4:9])+facet_wrap(~condition)+theme(strip.text=element_text(size=16),axis.text=element_text(size=12))
  ggsave(plt,filename = paste0(print.dir,'biplots/cond_biplot_',pc.use[1],ax.suf[1],'_',pc.use[2],ax.suf[2],'.pdf'),height=height,width=width)
}

make.mean.sigs <- function(data,cordir,num.pcs=10,levels=c('saline','coc.acute','coc.chron'),colors=c('blue','indianred2','brown'),
                           test.batch=T,test.max.p=T,max.p.bal=c(2,2),suffix='') {
  pc.scores=fread(paste0(cordir,'PC_allscores',suffix,'.csv'))
  pc.scores=pc.scores[,!"experiment",with=F]
  setkey(pc.scores)
  pc.scores=merge(pc.scores,unique(data[,.(cell.name,experiment)]),by='cell.name')
  pc.scores=merge(pc.scores,expt.cond,by='experiment')
  pc.scores[,cond.2:=condition];pc.scores[cond.2 %like% 'sal',cond.2:='saline'];
  pc.scores[,cond.2:=factor(cond.2,levels=levels)]

  pvals=data.table(pc.score=c(paste0('PC',1:num.pcs,'.score'),paste0('PC',1:num.pcs,'.pos'),paste0('PC',1:num.pcs,'.neg')))
  if(test.batch){
    pc.scores[condition %like% 'acute',batch:='acute'];pc.scores[condition %like% 'chron',batch:='chron']

  }
  setkey(pc.scores)
  for(pc in paste0('PC',1:num.pcs)){
    pc.scores[,mean.sig:=mean(get(paste0(pc,'.score'))),by=experiment]
    plt=plot.meansig(pc.scores,levels=levels,colors=colors,test.batch=test.batch,max.p.bal=max.p.bal)
    ggsave(plt[[1]],file=paste0(cordir,'coc_effects/',pc,'_meanSig_byExpt',suffix,'.pdf'),height=4,width=5)
    if(test.batch){
      pvals[pc.score==paste0(pc,'.score'),pbatch:=plt[[3]]]
    }
    pvals[pc.score==paste0(pc,'.score'),p.test:=plt[[2]]]

    pc.scores[,mean.sig:=mean(get(paste0(pc,'.pos'))),by=experiment]
    plt=plot.meansig(pc.scores,levels=levels,colors=colors,test.batch=test.batch,max.p.bal=max.p.bal)
    ggsave(plt[[1]],file=paste0(cordir,'coc_effects/',pc,'_meanPos_byExpt',suffix,'.pdf'),height=4,width=5)
    if(test.batch){
      pvals[pc.score==paste0(pc,'.pos'),pbatch:=plt[[3]]]
    }
    pvals[pc.score==paste0(pc,'.pos'),p.test:=plt[[2]]]

    pc.scores[,mean.sig:=mean(get(paste0(pc,'.neg'))),by=experiment]
    plt=plot.meansig(pc.scores,levels=levels,colors=colors,test.batch=test.batch,max.p.bal=max.p.bal)
    ggsave(plt[[1]],file=paste0(cordir,'coc_effects/',pc,'_meanNeg_byExpt',suffix,'.pdf'),height=4,width=5)
    if(test.batch){
      pvals[pc.score==paste0(pc,'.neg'),pbatch:=plt[[3]]]
    }
    pvals[pc.score==paste0(pc,'.neg'),p.test:=plt[[2]]]

  }
  write.csv(pvals,file=paste0(cordir,'coc_effects/pvals_sigs',suffix,'.csv'),quote=F,row.names = F)
}

plot.meansig <- function(pc.scores,levels=c('saline','coc.acute','coc.chron'),colors=c('blue','indianred2','brown'),
                         test.batch=T,test.max.p=T,max.p.bal=c(2,2)) {
  if(test.batch){
    stats.cells=unique(pc.scores[,.(mean.sig,experiment,condition,cond.2,batch)])
    test=aov(mean.sig ~ batch,data=stats.cells)
    p.batch=summary(test)[[1]][["Pr(>F)"]][1]

  }
  else{
    stats.cells=unique(pc.scores[,.(mean.sig,experiment,condition,cond.2)])
  }
  if(test.max.p){
    num.expts=length(unique(stats.cells[,experiment]))
    setorder(stats.cells,mean.sig)
    stats.cells[1:max.p.bal[1],cond.max:='top']
    stats.cells[(1+max.p.bal[1]):(max.p.bal[1]+max.p.bal[2]),cond.max:='bot']
    maxp.test=aov(mean.sig ~ cond.max,data=stats.cells)
    max.pval.1=summary(maxp.test)[[1]][["Pr(>F)"]][1]

    stats.cells[,cond.max:=NULL];setorder(stats.cells,mean.sig)
    stats.cells[1:max.p.bal[2],cond.max:='top']
    stats.cells[(1+max.p.bal[2]):(max.p.bal[1]+max.p.bal[2]),cond.max:='bot']
    maxp.test=aov(mean.sig ~ cond.max,data=stats.cells)
    max.pval.2=summary(maxp.test)[[1]][["Pr(>F)"]][1]

    max.pval=min(max.pval.1,max.pval.2)

  }
  max.sig=max(stats.cells[,mean.sig]);min.sig=min(stats.cells[,mean.sig]);
  if(length(unique(stats.cells[,cond.2]))==2){
    test=t.test(mean.sig ~ cond.2,data=stats.cells)
    pval=test$p.value
    test.type='Students t'
  }
  else{
    test=aov(mean.sig ~ cond.2,data=stats.cells)
    pval=summary(test)[[1]][["Pr(>F)"]][1]
    test.type='aov'
  }


  plt=ggplot(stats.cells,aes(cond.2,mean.sig,color=experiment,group=cond.2,fill=cond.2))+geom_boxplot(alpha=.5,outlier.shape = NA)+
    geom_point(size=3,position=position_jitter(width=.3))+
    scale_color_manual(values=expt2col)+theme_classic()+scale_fill_manual(values=colors)+
    annotate(geom='text',x = 2,y=max.sig,label=paste0(test.type,' p = ',format(pval,digits=2)))
  if(test.batch){
    plt=plt+annotate(geom='text',x = 2,y=.5*(max.sig+min.sig),label=paste0('batch p = ',format(p.batch,digits=2)))
  }
  if(test.max.p){
    plt=plt+annotate(geom='text',x = 2,y=min.sig,label=paste0('max p = ',format(max.pval,digits=2)))
  }

  # output
  if(test.batch){
    if(test.max.p){
      list(plt,pval,p.batch,max.pval)
    }
    else{
      list(plt,pval,p.batch)
    }
  }
  else{
    if(test.max.p){
      list(plt,pval,max.pval)
    }
    else{
      list(plt,pval)
    }
  }
}

make.mean.scores <- function(data,cordir,num.pcs=10,levels=c('saline','coc.acute','coc.chron'),colors=c('blue','indianred2','brown'),test.batch=T,suffix='') {
  pc.scores=fread(paste0(cordir,'PC_allscores',suffix,'.csv'))
  pc.scores=pc.scores[,c('cell.name',paste0('PC',1:num.pcs)),with=F]
  pc.scores=merge(pc.scores,unique(data[,.(cell.name,experiment)]),by='cell.name')
  pc.scores=merge(pc.scores,expt.cond,by='experiment')
  pc.scores[,cond.2:=condition];pc.scores[cond.2 %like% 'sal',cond.2:='saline'];
  pc.scores[,cond.2:=factor(cond.2,levels=levels)]

  pvals=data.table(pc.score=paste0('PC',1:num.pcs,''))
  if(test.batch){
    pc.scores[condition %like% 'acute',batch:='acute'];pc.scores[condition %like% 'chron',batch:='chron']
  }
  for(pc in paste0('PC',1:num.pcs)){
    pc.scores[,mean.sig:=mean(get(pc)),by=experiment]

    if(test.batch){
      stats.cells=unique(pc.scores[,.(mean.sig,experiment,condition,cond.2,batch)])
      test=aov(mean.sig ~ batch,data=stats.cells)
      p.batch=summary(test)[[1]][["Pr(>F)"]][1]
      pvals[pc.score==pc,pbatch:=p.batch]
    }
    else{
      stats.cells=unique(pc.scores[,.(mean.sig,experiment,condition,cond.2)])
    }
    max.sig=max(stats.cells[,mean.sig]);min.sig=min(stats.cells[,mean.sig]);

    if(length(unique(stats.cells[,cond.2]))==2){
      test=t.test(mean.sig ~ cond.2,data=stats.cells)
      pval=test$p.value
      test.type='Students t'
    }
    else{
      test=aov(mean.sig ~ cond.2,data=stats.cells)
      pval=summary(test)[[1]][["Pr(>F)"]][1]
      test.type='aov'
    }
    pvals[pc.score==pc,p.test:=pval]

    plt=ggplot(stats.cells,aes(cond.2,mean.sig,color=experiment,group=cond.2,fill=cond.2))+geom_boxplot(alpha=.5,outlier.shape = NA)+
      geom_point(size=3,position=position_jitter(width=.3))+
      scale_color_manual(values=expt2col)+theme_classic()+scale_fill_manual(values=colors)+
      annotate(geom='text',x = 2,y=max.sig,label=paste0(test.type,' p = ',format(pval,digits=2)))
    if(test.batch){
      plt=plt+annotate(geom='text',x = 2,y=min.sig,label=paste0('batch p = ',format(p.batch,digits=2)))
    }
    ggsave(plt,file=paste0(cordir,'coc_effects/',pc,'_meanPCscore_byExpt',suffix,'.pdf'),height=4,width=5)
  }
  write.csv(pvals,file=paste0(cordir,'coc_effects/pvals_scores.csv'),quote=F,row.names = F)
}
make.additive.sigs <- function(data,cordir,sigs=c('PC1.score','PC6.score'), func='add',levels=c('saline','coc.acute'),colors=c('blue','red'),test.batch = T){
  pc.scores=fread(paste0(cordir,'PC_allscores.csv'))
  pc.scores=merge(pc.scores,unique(data[,.(cell.name,experiment)]),by='cell.name')
  pc.scores=merge(pc.scores,expt.cond,by='experiment')

  if(test.batch){
    pc.scores[condition %like% 'acute',batch:='acute'];pc.scores[condition %like% 'chron',batch:='chron']
  }

  if(func=='add'){
    pc.scores[,mean.sig:=mean(get(sigs[1])+get(sigs[2])),by=experiment]
    pc.scores[,cond.2:=condition];pc.scores[cond.2 %like% 'sal',cond.2:='saline'];
    pc.scores[,cond.2:=factor(cond.2,levels=levels)]
    plt=plot.meansig(pc.scores,levels=levels,colors=colors,test.batch = test.batch)
    ggsave(plt,file=paste0(cordir,'coc_effects/',sigs[1],'+',sigs[2],'_mean_byExpt.pdf'),height=4,width=5)
  }
  if(func=='subtract'){
    pc.scores[,mean.sig:=mean(get(sigs[1])-get(sigs[2])),by=experiment]
    pc.scores[,cond.2:=condition];pc.scores[cond.2 %like% 'sal',cond.2:='saline'];
    pc.scores[,cond.2:=factor(cond.2,levels=levels)]
    plt=plot.meansig(pc.scores,levels=levels,colors=colors,test.batch = test.batch)
    ggsave(plt,file=paste0(cordir,'coc_effects/',sigs[1],'-',sigs[2],'_mean_byExpt.pdf'),height=4,width=5)
  }

}

cast.to.sigscores <- function(data,scoredir) {
  pc.scores <- fread(paste0(scoredir,'PCSigScores.csv'))
  cast.plot <- dcast.data.table(data,cell.name~gene,value.var='log10.cpm',fill = 0);cast.plot <- merge(cast.plot,pc.scores,by='cell.name')
  cast.plot
}

# Plot the means of top genes vs principal components to show what a discrete subtype looks like (step function)
# For custom x value, include a column with name matching the parameter x.custom
plot.gene.curves <- function(data,num.genes=40,pc.x='PC1.score',pc.genes='PC1',dir=cordir,num.to.avg=6,genes.custom='',x.custom='',height=28,width=34,method.clust='pearson',suffix='') {
  library(zoo)

  if(genes.custom==''){
    pc.loadings <- fread(paste0(dir,'PCloadings.csv'))
    if(length(pc.genes)==1){
      setorderv(pc.loadings,pc.genes); genes.pc1.pos <- pc.loadings[,gene][1:num.genes]
      setorderv(pc.loadings,pc.genes,order = -1); genes.pc1.neg <- pc.loadings[,gene][1:num.genes]
      genes.plot <- c(genes.pc1.pos,genes.pc1.neg)
    }
    else if(length(pc.genes)==2){
      setorderv(pc.loadings,pc.genes[1]); genes.pc1.pos <- pc.loadings[,gene][1:num.genes]
      setorderv(pc.loadings,pc.genes[1],order = -1); genes.pc1.neg <- pc.loadings[,gene][1:num.genes]
      setorderv(pc.loadings,pc.genes[2]); genes.pc2.pos <- pc.loadings[,gene][1:num.genes]
      setorderv(pc.loadings,pc.genes[2],order = -1); genes.pc2.neg <- pc.loadings[,gene][1:num.genes]
      genes.plot <- c(genes.pc1.pos,genes.pc1.neg,genes.pc2.pos,genes.pc2.neg)
      cast.plot=dcast.data.table(data[gene %in% genes.plot],cell.name~gene,value.var='log10.cpm',fill = 0)
    }
  }
  else{genes.plot=genes.custom}

  cast.plot=dcast.data.table(data[gene %in% genes.plot],cell.name~gene,value.var='log10.cpm',fill = 0)

  if(x.custom==''){
    pc.scores=fread(paste0(dir,'PC_allscores.csv'))
    x.vals=unique(pc.scores[,.(cell.name,get(pc.x))])
    setnames(x.vals,c('cell.name','x'))
  }
  else{
    x.vals=unique(data[,.(cell.name,x.custom)])
    setnames(x.vals,c('cell.name','x'))
    pc.x=x.custom
  }
  cast.plot=merge(cast.plot,x.vals,by='cell.name')
  cast.plot=data.frame(cast.plot[,2:ncol(cast.plot),with=F],row.names = cast.plot[,cell.name])
  cast.plot=cast.plot[order(cast.plot[,'x']),] # order data by x

  cast.means=as.data.table(apply(cast.plot,2,function (x) colMeans(matrix(x, nrow=num.to.avg))))
  data.rmean <- melt(cast.means,value.name = 'log10.cpm',id.vars = c('x'),variable.name = 'gene')
  data.rmean[,cell.name:=x]

  if(genes.custom[1]=='' | length(genes.custom) > 2){
    data.rmean[,gene:=factor(gene,levels=get.clustorder(data.rmean,method.cor=method.clust))]
  }



  plt=ggplot(data.rmean,aes_string('x','log10.cpm'))+geom_point(size=3,pch=1)+theme_classic()+facet_wrap(~gene)+
    theme(strip.text=element_text(size=20),axis.title=element_text(size=20))+geom_smooth(se=F,color='red')
  ggsave(plt,file=paste0(dir,'biplots/means_',pc.x,suffix,'.pdf'),height=height,width=width)

  data.rmean[,pc.order:=order(x),by=gene]
  plt=ggplot(data.rmean,aes_string('pc.order','log10.cpm'))+geom_point(size=3,pch=1)+theme_classic()+facet_wrap(~gene)+
    theme(strip.text=element_text(size=20),axis.title=element_text(size=20))+geom_smooth(se=F,color='red')
  ggsave(plt,file=paste0(dir,'biplots/means_',pc.x,suffix,'_order.pdf'),height=height,width=width)

  system(paste0('touch ',dir))
}

# function to get gene associated with a defined list of PC components
get.pc.genes <- function(pc.loadings, pcs.use, ng = 30){
  g.use = character(0)


  for(pc in pcs.use){


    if(nchar(pc) > 4) pc.short = substr(pc, 1, regexpr("[0-9]\\.", pc))
    else pc.short <- pc

    dc = F
    # sort genes in descending order if the "positive" component is requested, otherwise increasing order
    if(substr(pc, nchar(pc) - 2, nchar(pc)) == 'pos') {
      setorderv(pc.loadings, pc.short, -1)
      g.use = c(g.use, pc.loadings[, gene][1:ng])
    }
    else if(substr(pc, nchar(pc) - 2, nchar(pc)) == 'neg') {
      setorderv(pc.loadings, pc.short, 1)
      g.use = c(g.use, pc.loadings[, gene][1:ng])
    }
    else if(substr(pc, nchar(pc) - 2, nchar(pc)) == 'ore' | nchar(pc) < 5) {
      setorderv(pc.loadings, pc.short, -1)
      g.use = c(g.use, pc.loadings[, gene][1:ng])
      setorderv(pc.loadings, pc.short, 1)
      g.use = c(g.use, pc.loadings[, gene][1:ng])

    }
  }

  unique(g.use)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


theme_black = function(base_size = 12, base_family = "") {

  theme_grey(base_size = base_size, base_family = base_family) %+replace%

    theme(
      # Specify axis options
      axis.line = element_blank(),
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),
      axis.ticks = element_line(color = "white", size  =  0.2),
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),
      axis.ticks.length = unit(0.3, "lines"),
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),
      legend.key = element_rect(color = "white",  fill = "black"),
      legend.key.size = unit(1.2, "lines"),
      legend.key.height = NULL,
      legend.key.width = NULL,
      legend.text = element_text(size = base_size*0.8, color = "white"),
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),
      legend.position = "right",
      legend.text.align = NULL,
      legend.title.align = NULL,
      legend.direction = "vertical",
      legend.box = NULL,
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),
      panel.border = element_rect(fill = NA, color = "white"),
      panel.grid.major = element_line(color = "grey35"),
      panel.grid.minor = element_line(color = "grey20"),
      panel.margin = unit(0.5, "lines"),
      # Specify facetting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),
      strip.text.x = element_text(size = base_size*0.8, color = "white"),
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),
      plot.title = element_text(size = base_size*1.2, color = "white"),
      plot.margin = unit(rep(1, 4), "lines")

    )
}
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  frac.expr <- function (x, na.rm=FALSE) {
    if (na.rm) sum(x > 0)/length2(x,na.rm=T)
    else       sum(x > 0)/length2(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     frac = frac.expr    (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
## Summarizes data, handling within-subjects variables by removing inter-subject variability.
## It will still work if there are no within-S variables.
## Gives count, un-normed mean, normed mean (with same between-group mean),
##   standard deviation, standard error of the mean, and confidence interval.
## If there are within-subject variables, calculate adjusted values using method from Morey (2008).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   betweenvars: a vector containing names of columns that are between-subjects variables
##   withinvars: a vector containing names of columns that are within-subjects variables
##   idvar: the name of a column that identifies each subject (or matched subjects)
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySEwithin <- function(data=NULL, measurevar, betweenvars=NULL, withinvars=NULL,
                            idvar=NULL, na.rm=FALSE, conf.interval=.95, .drop=TRUE) {

  # Ensure that the betweenvars and withinvars are factors
  factorvars <- vapply(data[, c(betweenvars, withinvars), drop=FALSE],
                       FUN=is.factor, FUN.VALUE=logical(1))

  if (!all(factorvars)) {
    nonfactorvars <- names(factorvars)[!factorvars]
    message("Automatically converting the following non-factors to factors: ",
            paste(nonfactorvars, collapse = ", "))
    data[nonfactorvars] <- lapply(data[nonfactorvars], factor)
  }

  # Get the means from the un-normed data
  datac <- summarySE(data, measurevar, groupvars=c(betweenvars, withinvars),
                     na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Drop all the unused columns (these will be calculated with normed data)
  datac$sd <- NULL
  datac$se <- NULL
  datac$ci <- NULL

  # Norm each subject's data
  ndata <- normDataWithin(data, idvar, measurevar, betweenvars, na.rm, .drop=.drop)

  # This is the name of the new column
  measurevar_n <- paste(measurevar, "_norm", sep="")

  # Collapse the normed data - now we can treat between and within vars the same
  ndatac <- summarySE(ndata, measurevar_n, groupvars=c(betweenvars, withinvars),
                      na.rm=na.rm, conf.interval=conf.interval, .drop=.drop)

  # Apply correction from Morey (2008) to the standard error and confidence interval
  #  Get the product of the number of conditions of within-S variables
  nWithinGroups    <- prod(vapply(ndatac[,withinvars, drop=FALSE], FUN=nlevels,
                                  FUN.VALUE=numeric(1)))
  correctionFactor <- sqrt( nWithinGroups / (nWithinGroups-1) )

  # Apply the correction factor
  ndatac$sd <- ndatac$sd * correctionFactor
  ndatac$se <- ndatac$se * correctionFactor
  ndatac$ci <- ndatac$ci * correctionFactor

  # Combine the un-normed means with the normed results
  merge(datac, ndatac)
}
# pc.loadings <- fread('d2_corr4/PCloadings.csv')
# genes.plot <- unique(c(pc.loadings[order(PC9),gene][1:40],pc.loadings[order(-PC9),gene][1:40],
#                        pc.loadings[order(PC10),gene][1:40],pc.loadings[order(-PC10),gene][1:40]))
# pc.scores <- fread('d2_corr4/PCSigScores.csv')
# data.plot <- merge(data.d2[gene %in% genes.plot],pc.scores,by='cell.name')
# data.plot[,gene:=factor(gene,levels=get.clustorder(data.plot,method.cor='spearman'))]
# plt=ggplot(data.plot,aes(PC9.score,PC10.score,color=log10.cpm))+geom_point(size=2.2)+theme_classic()+
#   scale_color_gradientn(colours=brewer.pal(9,'YlOrBr')[4:9])+facet_wrap(~gene)+theme(strip.text=element_text(size=16))
# ggsave(plt,filename = 'd2_corr4/biplots/facet_wrap_PC9.score_PC10.score.pdf',height=28,width=34)


sem <- function(x) {sd(x)/sqrt(length(x))}

