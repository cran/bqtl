"summary.swap"<-
function(sc.obj, method=NULL, ncoef=length(sc.obj$alt.coef),nloc=sc.obj$nloc)
{
### summarize a swap.cycle object
  this.call <- sys.call()
  if (missing(method))
    method <-
      if(diff(dim(sc.obj$config))[1] < 0 )
        "F2"
      else
        "BC1" #includes RIL cases
  if(method=="F2") {
    if (is.null(ncoef)) ncoef <- ((max(sc.obj$config) - 1) %/% 2 + 1) * 2
    if (is.null(nloc)) nloc <- ncoef/ 2
  }
  else
    {
      if (is.null(ncoef)) {
        if (is.null(nloc))
          ncoef <- nloc <- max(sc.obj$config)
        else
          ncoef <- nloc
      }
      else
        nloc <- ncoef
    }
  coefs <- rep(0, ncoef)
  rep.post <- rep(sc.obj$post, rep(nrow(sc.obj$coefs), length(sc.obj$post
                                                              )))
  cp <- sc.obj$coefs * rep.post
  tap.coef <- tapply(cp, sc.obj$conf, sum)
  coefs[sort(unique(sc.obj$config))] <- if (method=="F2") tap.coef[-1] else tap.coef
  coefs <- coefs/sum(sc.obj$post)
  if (method=="F2") {
    locs <- apply((sc.obj$config - 1) %/% 2 + 1, 2:3, function(x) {unique(x[x > 0])})
  }
  else
    locs <- sc.obj$config
  diml <- dim(locs)
  if(length(diml) == 2)
    diml <- c(1, diml)
  loc.post <- rep(0, nloc)
  rep.marg <- rep(sc.obj$marg, rep(diml[1], diml[2] * diml[3]))
  loc.post[sort(unique(locs))] <-
    tapply(rep.marg, locs, sum)/sum(rep.marg)*diml[2]
  ## note assumption that second subscript tells size of model
  blk.ratio <- apply(sc.obj$cond/sc.obj$marg, 2, mean)
  res <- list(loc.posterior = loc.post, coefs = coefs,
              ratio = list(mean = mean(blk.ratio),
                var = var(blk.ratio)/length(blk.ratio)),
              setup=c(genes=diml[2],nreps=diml[3],nloc=nloc)
              ,call=this.call)
  class(res) <- "summary.swap"
  res
}
