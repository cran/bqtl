"summary.adj"<-
function(object, n.loc, coef.znames, mode.names = c("add", "dom"), imp.denom
	 = NULL, swap.obj = NULL,...)
{

### imp.denom is any non-null to switch on use of swap.obj$hk.exact
### this amounts to importance sampling from a swap.obj that differs
### in params from object to use enumerations (like one gene case)
### use imp.denom=1/elementwise.prior
###

  adj.unlist <-
    function(x,component="adj")
      sapply(x,"[[",component)
  adj.type <- if (is.null(imp.denom)) {
    if (is.null(swap.obj))
      "uniform"
    else
      "swap"
  }
  else {
    if (is.null(swap.obj))
      {
        imp.denom <- 1
        "uniform"
      }
    else
      "imp.swap"
  } 
  if(is.null(mode.names)) {
    coef.unames <- coef.znames
    n.mode <- 1
  }
  else {
    coef.unames <- outer(mode.names, coef.znames, paste, sep = ".")
    n.mode <- length(mode.names)
  }

  
  if(is.element(adj.type,c("swap","imp.swap"))) {
    config <- uniq.config(swap.obj)
    if(max(config$match) != length(object))
      stop("incompatible argument lengths:object and swap.obj"
           )
  }
  else
    config <- list(match=seq(along = object),
                   uniq=rbind(seq(along = object)))
  
  if (is.null(object[[1]]$reg.vec)) { # take locs from swap obj or default
    locs <- (config$uniq[,config$match] - 1)%/%n.mode +1
    coef.ind <- config$uniq[,config$match]
    coef.ind <- coef.ind[coef.ind>0]
  }
  else { # figure out locs from object[[i]]$reg.vec
    locs <-
      factor(unlist(lapply(object[config$match], function(x, nm, n.mode)
                           unique((match(x$reg.vec, nm) - 1) %/% n.mode) + 1,
                           nm = coef.unames, n.mode = n.mode)),
             1:length(coef.znames))
    coef.names <- unlist(lapply(object, "[[", "reg.vec")[config$match])
    coef.ind <- factor(coef.names, coef.unames)
  }
  coef.multiple <- unlist(lapply(object, function(x) length(x$parm) - 2))[config$match]
  
  coefs <- unlist(lapply(object, function(x) x$parm[ - c(1, length(x$parm))])[config$match])
  
  
  adj <- adj.unlist(object)[config$match]
  post <- adj*unlist(adj.unlist(object,component="hk.exact"))[config$match]
  switch(adj.type,
       uniform={
         post <- post/imp.denom
         adj.mean <- mean(adj/ (imp.denom/sum(imp.denom)) )
         adj.var <- 0
         hk.ratio.mean <- 1
         loc.post <- tapply(rep(post, rep(n.loc, length(post))),
                            locs, sum)
         loc.post[is.na(loc.post)] <- 0
         loc.post <- loc.post/sum(loc.post)
         coef.avg <- tapply(coefs*rep(post,coef.multiple), coef.ind, sum)/sum(post)
         coef.avg[is.na(coef.avg)] <- 0
         },
       swap={
         adj.mean <- mean(adj)
         adj.var <- var(c(apply(matrix(adj, nrow = n.loc), 2, mean)))/(
                                              length(adj)/n.loc)
         hk.ratio.mean <- 1
         loc.post <- tapply(rep(adj, rep(n.loc, length(adj))), locs, sum)
         loc.post[is.na(loc.post)] <- 0
         loc.post <- loc.post/sum(loc.post)
         
         coef.adj <- rep(adj, coef.multiple)
         coef.avg <- tapply(coef.adj * coefs, coef.ind, sum)/sum(adj)
         coef.avg[is.na(coef.avg)] <- 0
       },
       imp.swap={
         post <- post/swap.obj$hk.exact
         adj.mean <- mean(post)
         hk.ratio <- adj.unlist(object, comp="hk.exact")[config$match]/swap.obj$hk.exact
         hk.ratio.mean <- mean(hk.ratio)	
         adj.var <- var(c(apply(matrix(adj * hk.ratio, nrow = n.loc), 2, 
                                mean)))/(length(adj)/n.loc)
         loc.post <- tapply(rep(post, rep(n.loc, length(post))), locs, mean)
         loc.post[is.na(loc.post)] <- 1
         loc.post <- loc.post * swap.obj$alt.marg
         loc.post <- loc.post/sum(loc.post)
         coef.adj <- rep(post, coef.multiple)
         coef.avg <- tapply(coef.adj * coefs, coef.ind, sum)/sum(post)
         coef.avg[is.na(coef.avg)] <- 0
       }
       )
  list(adj = adj.mean, var = adj.var, coef = coef.avg, loc = loc.post, 
       hk.ratio.mean = hk.ratio.mean)
}
