library(Seurat)
library(mclust)
library(MASS)
library(parallel)
library(doMC)
library(dplyr)
library(FactoMineR)

# Notes 
# subsampling helps a lot (Arterial - cont w CV, SVc, SVv)
# Default, not sum-of-30, PC scores actually seem to "work" better - give more separation between cell types.
# Distributions don't seem more skewed necessarily
# negative correlation gene filter without subsampling does not seem to help
# what about downsampling cells to equal # per cell type? Seems to add noise, does not really help (or hurt)
# what about minimum cell ratio of 2?

# subsample + negative corr. + CPCA + sum-of-30: very large spread in PC1.neg and PC1.pos values even for very
# discrete cell types. Indication that CPCA is maximizing for variability *within* subtypes as well as between
# subtypes (this would be the expected behavior!)


#' Approximate PDF of transitioning/intermediate cells with a constant sampling rate.
#' for now, sig needs to be bounded above by < 0.17*(mu_a - mu_b) since the pdf gives an erroneous result above that. This bound on sigma (noise) is OK since it is difficult, and not very useful, to differentiate between two very wide Gaussians and a pure continuum/transition
#' @param x Numeric vector
#' @param sig Numeric, Width of the noise (standard deviation).     
#' @param mu_a Numeric, mean of beginning state
#' @param mu_b Numeric, mean of end state
#' 
#' @return A vector of probabilities of generating each input vector
#'
transit_pdf  <- function(x, sig=0.15, mu_a=-0.6, mu_b=0.6){    
  
  
  n=2.7 # these are tuned to make the approximate distribution (constant + Exp[(x-mu)^3] terms) close to the true distribution (convolution of step + Gaussians)
  w=0.93 # scaling the width of the Exp[(x-mu)^3] term relative to the sigma of the Gaussians for cell states
  
  exp_cubed <- function(y,mu,s){
    return(exp(-1*abs((y - mu)^3)/s^1.5))
  }
  exp_cubed2 <- function(y,mu,s){
    return(exp_cubed(y, mu, s) - exp_cubed(mu + 2*n*sig, mu, s))
  }
  exp_cubed3 <- function(y,mu,s){
    return(exp_cubed2(y, mu, s)*(exp_cubed2(mu, mu, s)^(-1)))
  }
  
  unscaled_density <- function(x, mu_a, mu_b, sig){
    return((x > mu_a - n*sig) * exp_cubed3(x, mu_a + n*sig, w*sig) * (x < mu_a + n*sig) + # noise
             (x >= mu_a + n*sig) * (x < mu_b - n*sig) + # constant
             (x >= mu_b - n*sig) * exp_cubed3(x, mu_b - n*sig, w*sig) * (x < mu_b + n*sig)) # noise
  }
  
  dx = abs((mu_a-n*sig - mu_b+n*sig)/1000)
  tscale = sum(dx*unscaled_density(seq(mu_a-n*sig, mu_b+n*sig, dx), mu_a=mu_a, mu_b=mu_b, sig=sig))^(-1)
  
  
  return(tscale * unscaled_density(x=x, mu_a=mu_a, mu_b=mu_b, sig=sig))
}

#' PDF of cells along a 1D axis of cell identity (A, B, and A-B intermediate).
#' The pdf of A-B intermediate cells are given by transit_pdf, which assumes 
#' constant sampling rate of intermediates across cell identity
#' The pdf of cells of type A or type B are given by Gaussians with equal widths
#' fa of the cells in state A, fab transitioning, and fb in state B
#' 
#' @return A vector of probabilities of generating each input vector
#'

cellid_pdf <- function(x, fa=0.5, fab=0.2, mu_a=-0.6, mu_b=0.6, sigma=0.1, max_sigma_scale=0.15){
  fb  = 1 - (fa + fab) # fraction of cells in type B
  
  if (mu_b < -1 | mu_a > 1
      | fa < 0 | fab < 0 | fb < 0 
      | sigma > max_sigma_scale*abs(mu_b-mu_a) # restrict noise term to force the fit of a "continuum" in 
      | fa > 1 | fb > 1 |
      (fa + fab + fb) > 1
      | sigma < 0.005*abs(mu_b-mu_a) # sigma needs to not be too small
      | mu_a > mu_b){ # constrain variables (for optimization)
    return(x*0)
  }
  else{
    # approximate normalization parameter for transitioning pdf
    return(fa * dnorm(x, mean=mu_a, sd=sigma) + 
             fb * dnorm(x, mean=mu_b, sd=sigma) + 
             fab * transit_pdf(x=x, sig=sigma, mu_a=mu_a, mu_b=mu_b))
  }
  
}

plot_cellid_fit <- function(data, params){
  
  print(params)
  xs = seq(-1, 1, 0.02);
  hist(data, col='red', breaks = 20, freq=F);
  lines(xs, cellid_pdf(xs, params[1], params[2], params[3], params[4], params[5]), 'l', xlab = 'x', ylab='P(x)')
}

fit_GaussianMM <- function(data){
  mcl.model <- Mclust(data, 2, warn = F, verbose = F)
  
  mu = as.numeric(mcl.model$parameters$mean)
  mu_a=mu[1]
  mu_b=mu[2]
  sigma = sqrt(mcl.model$parameters$variance$sigmasq)*.5
  f_a =sum(mcl.model$classification==1)/length(data)
  
  if(mu_a > mu_b){ # reverse guesses
    tmp  = mu_b
    mu_b = mu_a
    mu_a = tmp
    tmp  = f_a
    f_a = 1 - tmp
  }
  return(c(f_a, mu_a, mu_b, sigma))
}

calc_discrete_stat <- function(scores, n_rep=2){
  
  for(i in 1:n_rep){
    # estimate which cells are in the middle by fitting to a 2-Gaussian mixture model 
    res.2MM <- fit_GaussianMM(scores)
    scores_int = scores[(scores < res.2MM[3]) & (scores > res.2MM[2])]
    dists_int = sort(diff(sort(scores_int)), decreasing = T)
    
    # generate exponential distribution from the dists
    # use the largest 20 dists, except for the very largest dist (in case there is a very skewed 
    # cell fraction)
    dists_to_fit <- dists_int[2:min(length(dists_int), 30)]
    fit <- fitdistr(dists_to_fit, 'exponential', verbose=F)
    probs=log10(sort(dexp(dists_int, rate = fit$estimate)))
    
    test_stat=probs[1] - mean(probs[2:4])
    if(i==1){
      max_test=test_stat
    }
    else{
      max_test=max(max_test, test_stat)
    }
  }
  return(test_stat)
}

discrete_stat_bootstrap <- function(data, n_boot=100){
  
  stat_boot = numeric(0)
  for(i in 1:n_boot){
    data_boot = sample(data, replace = TRUE, size = length(data))
    stat_i = calc_discrete_stat(data_boot)
    stat_boot <- c(stat_boot, stat_i)
  }
  return(stat_boot)
}


# return(c(f_a, mu_a, mu_b, sigma))
nmfit_IntermediateScan <- function(data, n=5){
  
  # negative loglikelihood function to be minimized
  nll_cellid_IntScan <- function(params){
    return(-1*sum(log(cellid_pdf(x=data, fa=params[1], fab=params[2], mu_a=params[3], mu_b=params[4], sigma=params[5]))))
  }
  
  rebalance_fa <- function(fa, fab){
    # recalculate f_a so f_a/f_b stays the same and f_a+f_ab+f_b==1
    # assuming fab==0 initially (i.e., fa comes from fit_GaussianMM)
    fb = 1-fa
    alph = fa/(fa + fb)
    fa_new = (1 - fab) * fa/(1-fa) 
    return(fa_new)
  }
  
  i=1
  param_mins=matrix(nrow = 0, ncol = 5)
  ll_mins = numeric(0)
  for(fab in c(0, seq(0.05, .98, length.out = n))){
    mm_params = fit_GaussianMM(data)
    fa_init <- (1 - fab) * mm_params[1]
    
    # don't let sigma guess be too large
    sigma_init = pmax(0, pmin(0.15*abs(mm_params[2]-mm_params[3]), mm_params[4])) 
    x0 = c(fa_init, fab, mm_params[2:3], sigma_init)
    #print("Initial Guess")
    #print(x0)
    optim.res=stats::optim(par=x0, fn = nll_cellid_IntScan)
    xopt = optim.res[1]$par
    nll  = optim.res[2]$value
    #print(xopt)
    #print(nll)
    #print(xopt)
    
    ll_mins    = c(ll_mins, nll)
    param_mins = rbind(param_mins, xopt)
    
  }
  return(list(ll_mins, param_mins))
}

fit_cellIdentityDistr <- function(scores, n=5){
  
  fit.res = nmfit_IntermediateScan(scores, n=n)
  ll = fit.res[[1]]
  params = fit.res[[2]]
  
  argmin <- function(x) return(which(x==min(x)))
  
  return(params[argmin(ll),])
}



#' Perform pairwise continuity analysis on grouping of cells Possible: sim genes
#' by correlation? pairwiseDiscretenessTest_v2: no longer calculating pvalues
#' for each gene (takes forever) 4/8/18: Implementing new cell scoring measure:
#' use the most discrete results from signatures based on DE genes (as in V2)
#' and CPCA with sum-of-30 modified PC scores.
RunPDT <- function(object, subsample.genes=F,
                                          num.sig.genes=30, max.cell.ratio=5, genes.use=NULL, dir=getwd(),
                                          n.cores=NULL) {
  
  
  
  if(!file.exists(dir)) system(paste("mkdir", dir))
  if(is.null(genes.use)) genes.use <- rownames(object@data)
  if(is.null(n.cores)) n.cores <- detectCores()
  if(is.na(n.cores)) n.cores <- 1
  #if(!exists(dir))  system(paste0('mkdir ', dir))
  cellids <- as.character(object@ident)
  ids.iterator <- make.IDpair.iterator(cellids)
  print(ids.iterator)
  # res <- foreach(i = 1:nrow(ids.iterator), .combine = "c") %do% {
  res=list()
  registerDoMC(n.cores)
  # calculate pairwise statistics
  pair.res <- foreach(i = 1:nrow(ids.iterator), .combine="c") %dopar% {
    id1=ids.iterator[i, 1]
    id2=ids.iterator[i, 2]
    
    numcells.use <- max.cell.ratio*min(sum(cellids==id1), sum(cellids==id2))
    
    pairObject <- SubsetData(object, ident.use = c(id1, id2), max.cells.per.ident = numcells.use,do.center = F,subset.raw = T)
    
    
    ident.rest <- scoreCellsByIdentitySignature(object=pairObject, id1 = id1, id2 = id2, num.sig.genes = num.sig.genes,
                                               genes.use=genes.use, dir=dir)
    
    
    # if there are not enough differentially expr genes, scoreCellsByIdentitySignature returns NULL
    # the cell types are then considered to be completely continuous as they cannot be separated.
    # Future versions will report the # of differentially expressed genes as a metric of the distance
    # between cell types. 
    if(is.null(ident.res)){
      dstats=rep(0, 50)
      cellid.params=c(0,1,-.5,.5,.1)
      ident.signature=data.frame(sig.1=rep(0,length(pairObject@cell.names)),
                                 sig.2=rep(0,length(pairObject@cell.names)),
                                 score=rep(0,length(pairObject@cell.names)))
    }
    
    # if there are sufficient diff expr genes, compare the diff expr scores to
    # scores generated by PCA and use the one that is more discrete. 
    else{
      
      pc.res <- scoreCellsByBestPC(object = pairObject, id1 = id1, id2 = id2, n.pcs = 3, genes.use = genes.use)
      
      
      
      # calculate the discrete statistics and fit the 5-parameter model of cell
      # identity (2-Gaussian + intermediates)
      ident.signature <- ident.res[[1]]
      
      dstats <- discrete_stat_bootstrap(ident.signature$score, n_boot = 50)
      cellid.params <- fit_cellIdentityDistr(ident.signature$score, n = 5)
    }
    

    
    # genes used to calculate signature
    genes.1 <- ident.res[[2]]
    genes.2 <- ident.res[[3]]
    
    list(c(id1, id2), dstats, cellid.params, ident.signature)#, genes.1, genes.2)
    
  }
  
  # summarize for output
  discrete_thresh <- -4
  id.pairs         <- as.data.frame(ids.iterator)
  colnames(id.pairs) <- c("id1","id2")
  signature.scores <- list()
  # signature.genes.1 <- list()
  # signature.genes.2 <- list()
  for(i in 0:(nrow(id.pairs)-1)){
    discrete_bootstrap <- pair.res[[4*i + 2]]
    print(head(discrete_bootstrap))
    id.pairs$discrete.median[i+1] <- median(discrete_bootstrap)
    id.pairs$discrete.conf[i+1] <- sum(discrete_bootstrap < discrete_thresh)/length(discrete_bootstrap)
    
    id1=id.pairs$id1[i+1]
    id2=id.pairs$id2[i+1]
    ss.i <- pair.res[[4*i + 4]]
    colnames(ss.i) <- c(paste(id1, id2, "sig.1", sep="_"),
                        paste(id1, id2, "sig.2", sep="_"),
                        paste(id1, id2, "score", sep="_"))
    signature.scores <- c(signature.scores, ss.i)
    
    transitionFit_params <- pair.res[[4*i + 3]]
    id.pairs$f.transition[i+1] <- transitionFit_params[2]
    id.pairs$N.transition[i+1] <- transitionFit_params[2] * nrow(ss.i)
    
  }
  
  id.pairs$continuous.conf <- 1 - id.pairs$discrete.conf
  pair.res <- list(pairwise.statistics=id.pairs, 
                   signature.scores=signature.scores, 
                   cell.names=object@cell.names)
  save(pair.res, file=file.path(dir,"pairwise_results.Rdata"))
  return(pair.res)
  
}

#' Make a matrix of all pairwise combinations of ids
make.IDpair.iterator <- function(cellids){
  ids <- unique(as.character(cellids))
  return(t(combn(ids, m=2)))
}


# scoreCellsByPCA <- function(object, id1, id2, num.sig.genes=30, method = 'pct_or_diff', 
#                             genes.use=NULL, dir=getwd()){
#   
#   if(is.null(genes.use)) genes.use <- rownames(object@data)
#   
# }

#' Score cells by the sum of two sets of genes: genes.1 enriched in id1, and genes.2 enriched in id2
#' #' @param object Seurat object containing cells with id1 and id2 in the ident field
#' pairwiseDiscretenessTest_v2: no longer calculating pvalues for each gene (takes forever)
#' Change 4/6/18: Fixed a fault that gave the incorrect genes for calculating the id1 signature
scoreCellsByIdentitySignature <- function(object, id1, id2, num.sig.genes=30, method = 'pct_or_diff', genes.use=NULL){
  
  if(is.null(genes.use)) genes.use <- rownames(object@data)
  
  # subset genes
  
  # Get top genes by fold-change and differential expression
  fc <- sortGenes.logFC(object = object,ident.1 = id1, ident.2 = id2, genes.use = genes.use)
  pct <- sortGenes.logFC(object = object,ident.1 = id1, ident.2 = id2, genes.use = genes.use)
  
  genes.fc.1 <- names(fc)[fc>0]
  genes.pct.1 <- names(pct)[pct>0]
  
  # reverse the order of the differential expression values to get genes higher in cell type 2
  fc.rev <- fc[length(fc):1]
  pct.rev <- pct[length(pct):1]
  genes.fc.2 <- names(fc.rev)[fc.rev<0]
  genes.pct.2 <-  names(pct.rev)[pct.rev<0]
  
  
  genes.1 <- unique(c(genes.fc.1[1:min(num.sig.genes, length(genes.fc.1))],
                      genes.pct.1[1:min(num.sig.genes, length(genes.pct.1))]))
  
  genes.2 <- unique(c(genes.fc.2[1:min(num.sig.genes, length(genes.fc.2))],
                      genes.pct.2[1:min(num.sig.genes, length(genes.pct.2))]))
  
  
  # print(genes.1)
  # print(genes.2)
  # if there are not enough differentially expressed genes between subtypes, 
  if(min(length(genes.1), length(genes.2)) < 5){
    return(NULL)
  }
  else{
    sig.scores <- calc.signature(SubsetData(object, ident.use = c(id1, id2)), genes.1, genes.2)
    return(list(sig.scores, genes.1, genes.2))
  }
  
}



#' Do PCA on a pair of cell groups (id1 and id2). Choose the best dimension as
#' the one that most discretely separates cells into groups. 
scoreCellsByBestPC <- function(object, id1, id2, n.pcs = 3, pca.type = 'CPCA', genes.use=NULL){
  
  if(is.null(genes.use)) genes.use <- rownames(object@data)
  
  scorelist=list()
  pca <- do.pca(object, ncp = n.pcs, genes.use = genes.use)
  disc.stats <- c()
  for(i in 1:n.pcs){
    scores <- calc.sumof30.scores(object = object, loadings = pca[[3]], dim = i) 
    scorelist <- c(scorelist, scores)
    disc.stats <- c(disc.stats, calc_discrete_stat(scores=scores[, 3]))
  }

  print(disc.stats)
  dim.best <- which.min(disc.stats)
  scores.best <- calc.sumof30.scores(object = object, loadings = pca[[3]], dim = dim.best) 
  genes.best <- get.top30.genes(loadings=pca[[3]], dim=dim.best)
  
  return(list(scores.best, genes.best[[1]], genes.best[[2]]))
  
}

do.pca <- function(object,ncp=3, genes.use=NULL){
  
  if(is.null(genes.use)) {
    object <- FindVariableGenes(object,do.plot = F, display.progress = F)
    genes.use <- object@var.genes

  }
  cast.counts <- t(as.matrix(object@data[genes.use, ]))
  
  res <- PCA(cast.counts, ncp=ncp, graph=F,scale.unit = F)
  coord <- as.data.frame(res$ind$coord[, 1:ncp])
  colnames(coord) <- paste0("PC", 1:ncp)
  coord$cell.name <- rownames(coord)
  
  loadings <- as.data.frame(res$var$coord[, 1:ncp])
  colnames(loadings) <- paste0("PC", 1:ncp)
  loadings$gene <- rownames(loadings)
  
  return(list(res,coord,loadings))
}

calc.sumof30.scores <- function(object, loadings, dim, num.sig.genes=30){

  top30 <- get.top30.genes(loadings = loadings, dim = dim)
  
  scores <-get.sig(object = object, genes.pos = top30[[1]], genes.neg = top30[[2]])
  scores <- scores[,c('pos','neg','score')]
  colnames(scores) <- c(paste0('PC',dim,'.pos'),paste0('PC',dim,'.neg'), paste0('PC',dim,'.score'))
  return(scores)
}

get.top30.genes <- function(loadings, dim, num.sig.genes=30){
  genes.pos <- loadings$gene[order(loadings[, paste0('PC',dim)], decreasing = T)][1:num.sig.genes]
  genes.neg <- loadings$gene[order(loadings[, paste0('PC',dim)])][1:num.sig.genes]
  return(list(as.character(genes.pos), as.character(genes.neg)))
}

get.sig <- function(object,genes.pos,genes.neg){
  # generate signature value of top PC +/- correlated genes
  genes.use <- c(genes.pos,genes.neg)
  
  # print(genes.use)
  data <- t(as.matrix(object@data[genes.use, ]))
  data.scale = data / apply(data, 2, max)
  

  cast.sig <- data.frame(neg = rowSums(data.scale[, genes.neg]),
                         pos = rowSums(data.scale[, genes.pos]))
  rownames(cast.sig) <- rownames(data.scale)
  cast.sig$pos <- cast.sig$pos / max(cast.sig$pos)
  cast.sig$neg <- cast.sig$neg / max(cast.sig$neg)
  cast.sig$score <- cast.sig$pos - cast.sig$neg
  cast.sig
  
}


# Change 4/6/18: output is ordered in decreasing differential expresssion values
sortGenes.logFC <- function(object, ident.1, ident.2, genes.use=NULL, thresh=0.2){
  if(is.null(genes.use)) genes.use = rownames(object@data)
  else genes.use = genes.use[genes.use %in% rownames(object@data)]
  
  data.1 = t(as.matrix(object@data[genes.use, object@ident==ident.1]))
  data.2 = t(as.matrix(object@data[genes.use, object@ident==ident.2]))
  mean.1 = apply(data.1, 2, mean)
  mean.2 = apply(data.2, 2, mean)
  
  diff <- mean.1 - mean.2
  diff <- diff[abs(diff) > thresh]
  return(diff[order(diff, decreasing = T)])
  
}

sortGenes.pctDiff <- function(object, ident.1, ident.2, genes.use=NULL, thresh=0.05){
  if(is.null(genes.use)) genes.use = rownames(object@data)
  else genes.use = genes.use[genes.use %in% rownames(object@data)]
  
  data.1 = t(as.matrix(object@data[genes.use, object@ident==ident.1]))
  data.2 = t(as.matrix(object@data[genes.use, object@ident==ident.2]))
  pct.1 = apply(data.1 > 0, 2, sum) / nrow(data.1)
  pct.2 = apply(data.2 > 0, 2, sum)/ nrow(data.2)
  
  diff <- pct.1 - pct.2
  diff <- diff[abs(diff) > thresh]
  return(diff[order(diff, decreasing = T)])
  
}
#' Generates two scores ('pos' and 'neg') for cells based on two sets of genes
calc.signature <- function(object,genes.sig.1,genes.sig.2, method.scale='max'){
  
  genes.use <- c(genes.sig.1,genes.sig.2)
  data <- t(as.matrix(object@data[genes.use, ]))
  data.scale = data / apply(data, 2, max)
  
  cast.sig <- data.frame(sig.2 = rowSums(data.scale[, genes.sig.2]),
                         sig.1 = rowSums(data.scale[, genes.sig.1]))
  rownames(cast.sig) <- rownames(data.scale)
  cast.sig$sig.1 <- cast.sig$sig.1 / max(cast.sig$sig.1)
  cast.sig$sig.2 <- cast.sig$sig.2 / max(cast.sig$sig.2)
  cast.sig$score <- cast.sig$sig.1 - cast.sig$sig.2
  cast.sig <- cast.sig[, c("sig.1", "sig.2", "score")]
  
  return(cast.sig)
  
}

gene.feature.correlations <- function(object, feature.use = 'ident'){
  if(feature.use=='ident'){
    return(cor(object@ident, object@scale.data))
  }
}

# connectogram: overlay f_ab graph onto tSNE plot
PlotConnectogram <- function(object, pair.res, colorvec=NULL, median.disc.cutoff = -6,
                             disc.conf.cutoff = 0.6,  edge.width.type = 'fractionOfSmallest.intermediate'){
  
  require(igraph)
  if(is.null(colorvec)) colorvec <- gg_color_hue(length(unique(object@ident)))[as.factor(object@ident)]
  
  pairwise.stats <- pair.res$pairwise.statistics
  ids = unique(as.character(object@ident))
  nodes = as.data.frame(t(sapply(ids, function(x) apply(GetCellEmbeddings(object, 'tsne', cells.use=object@cell.names[object@ident==x]), 2, median))))
  nodes$ident <- rownames(nodes)                                   
  nodes <- merge(nodes, unique(data.frame(ident=object@ident, color=colorvec)), by = 'ident')
  nc.id <- as.data.frame(table(object@ident))
  colnames(nc.id) <- c('ident','numcells')
  nodes <- merge(nodes, nc.id, by='ident')
  
  
  edges = pairwise.stats[, c("id1","id2","discrete.median","discrete.conf","f.transition", "N.transition")]
  edges = merge(edges, data.frame(id1 = nodes$ident, x0 = nodes$tSNE_1, y0 = nodes$tSNE_2), by='id1')
  edges = merge(edges, data.frame(id2 = nodes$ident, x1 = nodes$tSNE_1, y1 = nodes$tSNE_2), by='id2')
  edges = edges[edges$f.transition > 0, ]
  edges = edges[edges$discrete.median > median.disc.cutoff, ]
  edges = edges[edges$discrete.conf < disc.conf.cutoff, ]
  
  edges$alph <- edges$discrete.median
  edges$alph[edges$discrete.median > -1] <- -1
  edges$alph[edges$discrete.median < median.disc.cutoff] <- median.disc.cutoff
  edges$alph <- (edges$alph - min(edges$alph))/(max(edges$alph) - min(edges$alph))
  
  
  # edges$sim <- 1-edges$discrete.conf
  
  # edges$sim <- edges$discrete.median
  # edges$sim[edges$discrete.median > -1] <- -1
  # edges$sim[edges$discrete.median < median.disc.cutoff] <- median.disc.cutoff
  # edges$sim <- (edges$sim - min(edges$sim))/(max(edges$sim) - min(edges$sim))
  # edges$sim[edges$discrete.conf > disc.conf.cutoff] <- 0
  # 
  # 
  # edges.2 <- edges
  # tmp <- edges.2$id2
  # edges.2$id2 <- edges.2$id1
  # edges.2$id1 <- tmp
  # edges.3 <- rbind(edges[,c("id1","id2","sim")], edges.2[,c("id1","id2","sim")])
  # edge.mat <- dcast(edges.3, id1 ~ id2, value.var = 'sim', fill = 0)
  # rownames(edge.mat) <- edge.mat$id1
  # sim.mat <- as.matrix(edge.mat[, 2:ncol(edge.mat)])
  # sim.mat <- sim.mat[, rownames(sim.mat)]
  # 
  # g = graph_from_adjacency_matrix(sim.mat, mode='undirected', weighted = T,diag = F)
  # 
  # vtx.col <- unique(data.frame(ident = as.character(object@ident),
  #                              color = as.character(colorvec)))
  # rownames(vtx.col) <- vtx.col$ident
  # 
  # vtx.col <- vtx.col[V(g)$name, c("color")]
  # 
  # plot.igraph(g, vertex.size=6,
  #             vertex.color=vtx.col,
  #             edge.color=alpha('black',alpha=E(graph = g)$sim),
  #             vertex.label=NA,
  #             edge.width = 2*E(graph = g)$sim)
  # 
  # return(g)
  
  
  
  if(edge.width.type == 'fraction.intermediate'){
    edges$width <- (10*(edges$f.transition))
  }
  else if(edge.width.type == 'number.intermediate'){
    edges$width <- .1*(edges$N.transition)
  }
  else if(edge.width.type == 'fractionOfSmallest.intermediate'){
    nc.id1 <- nc.id
    colnames(nc.id1) <- c("id1", "numcells.id1")
    nc.id2 <- nc.id
    colnames(nc.id2) <- c("id2", "numcells.id2")
    edges <- merge(edges, nc.id1, by="id1")
    edges <- merge(edges, nc.id2, by="id2")
    edges$minNcells <- min(edges$numcells.id1, edges$numcells.id2)
    # print(head(edges))
    
    edges$width <- sqrt((edges$N.transition / edges$minNcells))*5
  }
  else{
    edges$width <- 3
  }
  plot(object@dr$tsne@cell.embeddings, col = alpha(colorvec,.2), pch=20, asp=1) +
    with(edges, segments(x0,y0,x1,y1, lwd = width, col = alpha('black', edges$alph))) +
    points(nodes$tSNE_1, nodes$tSNE_2, col=as.character(nodes$color), cex=sqrt(nodes$numcells)*.35, pch=19)
  points(nodes$tSNE_1, nodes$tSNE_2, col="black", cex=sqrt(nodes$numcells)*.35, pch=1)
  
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
