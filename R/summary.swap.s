"summary.swap"<-
function(object, method=NULL, ncoef=length(object$alt.coef),nloc=object$nloc,...)
{
### summarize a swap.cycle object
  this.call <- sys.call()
  if (exists("unique.default")) unique <- unique.default
  if (missing(method))
    method <-
      if(diff(dim(object$config))[1] < 0 )
        "F2"
      else
        "BC1" #includes RIL cases
  if(method=="F2") {
    if (is.null(ncoef)) ncoef <- ((max(object$config) - 1) %/% 2 + 1) * 2
    if (is.null(nloc)) nloc <- ncoef/ 2
  }
  else
    {
      if (is.null(ncoef)) {
        if (is.null(nloc))
          ncoef <- nloc <- max(object$config)
        else
          ncoef <- nloc
      }
      else
        nloc <- ncoef
    }
  coefs <- rep(0, ncoef)
  rep.post <- rep(object$post, rep(nrow(object$coefs), length(object$post
                                                              )))
  cp <- object$coefs * rep.post
  tap.coef <- tapply(cp, object$conf, sum)
  coefs[sort(unique(c( object$config )))] <- if (method=="F2") tap.coef[-1] else tap.coef
  coefs <- coefs/sum(object$post)
  if (method=="F2") {
    locs <- apply((object$config - 1) %/% 2 + 1, 2:3, function(x) {unique(c( x[x > 0]) )})
  }
  else
    locs <- object$config
  diml <- dim(locs)
  if(length(diml) == 2)
    diml <- c(1, diml)
  loc.post <- rep(0, nloc)
  rep.marg <- rep(object$marg, rep(diml[1], diml[2] * diml[3]))
  loc.post[sort(unique(c( locs )))] <-
    tapply(rep.marg, locs, sum)/sum(rep.marg)*diml[2]
  ## note assumption that second subscript tells size of model
  blk.ratio <- apply(object$cond/object$marg, 2, mean)
  res <- list(loc.posterior = loc.post, coefs = coefs,
              ratio = list(mean = mean(blk.ratio),
                var = var(blk.ratio)/length(blk.ratio)),
              setup=c(genes=diml[2],nreps=diml[3],nloc=nloc)
              ,call=this.call)
  class(res) <- "summary.swap"
  res
}
