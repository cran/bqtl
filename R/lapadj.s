"lapadj"<-
function(reg.formula, ana.obj,  
         rparm = NULL,  tol = 9.900000000000001e-09,
	return.hess = FALSE, mode.names = NULL, mode.mat = NULL,
          maxit = 100, nem = 1,setup.only=FALSE,subset=NULL,casewt=NULL,
         start.parm=NULL, ...)
{
  this.call <- match.call(expand.dots=TRUE)
  n.data<-nrow(ana.obj$data)
  if (class(ana.obj)!="analysis.object") stop("ana.obj must be an analysis.object") 

  if (is.character(ana.obj$method))
    crstype <- if (is.element(ana.obj$method,c("BC1","RI.self","RI.sib"))) 1 else 2

  if(is.null(mode.mat))
      mode.mat <- ana.obj$mode.mat
  if(is.null(mode.names))
      mode.names <- dimnames(mode.mat)[[2]]
  
  iter <- maxit
  regr.names <- ana.obj$reg.names
  form.terms <- terms(reg.formula[-2],specials="covar")
  form.factors <- attr(form.terms,"factors")
  if (length(form.factors) == 0) { # intercept only, assume rparm[1]==0
      
      mf.expr <- expression(model.frame(reg.formula,ana.obj$data,
                        na.action=na.omit,subset=subset))[[1]]
      mf.expr$subset<- subset
      mf <- eval(mf.expr)
      y <-
        model.extract(   mf,                     "response")
      n  <- length(y)
      mean.y <- mean(y)
      var.y <- var(y)*(n-1)/n
      logpost <- sum(log(dnorm(y, mean.y , sqrt(var.y))))
      post <-
          exp((-(n - 1)/2*log(2*pi)-log(n)/2 -
               (n-1)/2*log(var.y*(n - 1)/2) +
               lgamma((n - 1)/2))-log(2))
      parm <- c(mean.y,log(var.y)/2)
      return(list(logpost=logpost,posterior=post,hk.exact=post,hk.approx=post,
                  parm=parm,reg.vec=character(0),call=this.call))
      ## return the null posterior
  }
  e.marker <- seq(nrow(form.factors)) %in% attr(form.terms,"specials")$covar
  actual.covar <- ifelse(e.marker, FALSE,
                         !dimnames(form.factors)[[1]] %in% ana.obj$reg.names)
  covar.ind <- ( e.marker | actual.covar )
  covar.terms <- covar.ind%*%form.factors>0
  plain.terms <- (1-covar.ind)%*%form.factors>0
  both.terms <- covar.terms & plain.terms
  if (any(both.terms)) {
      covar.terms <- ( !both.terms ) & covar.terms
      plain.terms <- ( !both.terms ) & plain.terms
      both.index <- seq(along=both.terms)[both.terms]
      both.covar <- form.factors[covar.ind,both.terms] %equiv% form.factors[covar.ind,covar.terms]
      if (!apply(both.covar,1,any)) stop("can't find covariate main effect for interaction")
      both.plain <- form.factors[!covar.ind,both.terms] %equiv% form.factors[!covar.ind,plain.terms]
      if (!apply(both.plain,1,any)) stop("can't find locus main effect for interaction")
  }
  
  if (any(covar.ind)) {
      lhs <- paste(deparse(reg.formula[[2]]),"~")
      covar.formula <-
          eval(parse(text=paste(c(lhs,dimnames(form.factors)[[2]][covar.terms]),
                     collapse="+")))
      covar.terms <- terms(covar.formula)
      cfe <- expression(model.frame(covar.terms,ana.obj$data,
                                 na.action=na.omit,subset=subset))[[1]]
      cfe$subset <- subset
      covar.frame <- eval(cfe)
      subset <- row.names(ana.obj$data)%in%row.names(covar.frame)
      covar.matrix <- model.matrix(covar.terms,covar.frame)
      y <- model.extract(covar.frame,"response")
      covar.assign <-
          if (exists("is.R") && is.R() )
              attr(covar.matrix,"assign")[-1]
          else
              rep(seq(along=attr(covar.matrix,"assign"))-1,sapply(attr(covar.matrix,"assign"),length))[-1] 

      covar.reps <- table(factor(covar.assign,sort(unique(covar.assign))))
      if (any(both.terms)) {
          pttz <- apply(both.covar,1,function(x) seq(x)[x])-1
          pttx <- apply(both.plain,1,function(x) seq(x)[x])-1
          nt2uz <- length(pttx)
      }
      else {
          pttz <- pttx <- nt2uz <- 0
      }
      ncz <- ncol(covar.matrix)-1
      nz2uz <- ncz
      ptz <- seq(along=covar.assign)-1 
  }
  else {
      y <- eval(reg.formula[[2]],ana.obj$data)
      
      if (is.null(subset))
        subset <- !is.na(y)
      if (length(y)!=nrow(ana.obj$loc.right) )
        stop("length of dependent variable does not match other args")
      y <- y[subset]
      covar.matrix <- NULL
      ptz <- pttx <- pttz <- nt2uz <- nz2uz <- ncz <- 0
  }
  
  if (any(plain.terms)) {
      plain.formula <-
          eval(parse(text=paste(c("~",dimnames(form.factors)[[2]][plain.terms]),collapse="+")))
      reg.vec <- if (exists("is.R") && is.R() )
          as.character(attr(form.terms,"variables")[-1])[plain.terms]
      else
          as.character(attr(form.terms,"variables"))[plain.terms]
      
      which.vars <- sort(match(reg.vec, regr.names))
      
      if (length(which.vars)==0)
          stop("regressors not recognized - check spelling")
      ncx <- length(which.vars)
      cols.used <- col(regr.names)[which.vars]
      z.used <- unique(cols.used)
      nloc <- length(z.used)
      nrx <- if (crstype==1)
          2^nloc
      else
          3^nloc
      
      if (crstype==1) {
          if (nloc==1) reg.frame <- as.data.frame(mode.mat)
          else
              reg.frame <- as.data.frame(do.call("expand.grid", rep(list(mode.mat), nloc)))
          marker.distances <- switch(ana.obj$method,
                                     "RI.self"={ana.obj$map.frame$lambda/(2-ana.obj$map.frame$lambda)},
                                     "RI.sib"={ana.obj$map.frame$lambda/(4-3*ana.obj$map.frame$lambda)},
                                     ana.obj$map.frame$lambda)
      }
      else
      {
          marker.distances <- ana.obj$map.frame$lambda
          modes.used <- row(regr.names)[which.vars]
          col.seq <- match(cols.used,z.used)
          loc.indx <- as.matrix(do.call("expand.grid",rep(list(1:3),nloc)))
          reg.frame <-
              as.data.frame(matrix(mode.mat[cbind(c(loc.indx[,col.seq]),
                                                  rep(modes.used,rep(nrx,ncx)))],
                                   nc=ncx))
      }
      names(reg.frame) <- regr.names[which.vars]
      reg.proto <- model.matrix(plain.formula,reg.frame)[,-1,drop=FALSE]
      ncx <- ncol(reg.proto)
      ptx <- seq(from=0,by=1,length=ncx)
  } # if (any(plain.terms)) {
  else
  {
      nrx <- 1
      ncx <- nloc <- reg.proto <- ptx <-  z.used <- modes.used <- 0
  }

### set up nparm
  nx2uz <- ncx
  nreg <- 1 + nx2uz + nz2uz + nt2uz
  np <- nreg+1
  
  n <- length(y)
  if (is.null(start.parm)) { # initialize coefs and ln.sigma?
    ninit <- 0
    coefs <- double(nreg)
    ln.sigma <- double(1)
  }
  else {
    ninit <- 1
    if (length(start.parm)!=np) stop("length(start.parm) is incorrect")
    coefs <- as.double(start.parm[-np])
    ln.sigma <- as.double(start.parm[np])
  }
  nparm <- c(n,nrx,ncx,nloc,ncz,
             nx2uz,nz2uz,nt2uz,
             nreg,np,ptx,ptz,pttx,pttz,ninit)
  reg.vec <- c(dimnames(reg.proto)[[2]],dimnames(covar.matrix)[[2]][-1],dimnames(form.factors)[[2]][both.terms])
  
### rparm for intercept is assumed to be in the first element, named "intercept",
###  and names(rparm) must be matched, except when length matches number of variables
###  in regr.names
### interactions can ONLY be handled by scalar rparm unless
### names(rparm)%in%dimnames(reg.proto)[[2]]
###
  rparm.named <- !is.null(names(rparm))
  
  if(is.null(rparm)) {
      rparm <- rep(0, nreg)
  }
  else {
      if (rparm.named)
          rparm <- c(rparm["intercept"],rparm[reg.vec])
      else
      {
          len.choices <- c(1,length(regr.names)+1,nreg,2)
          
          rparm <-
              switch(match(length(rparm),len.choices,5),
                     c(0,rep(rparm,nreg-1)),
                 {
                     names(rparm) <- c("intercept",regr.names)
                     c(rparm["intercept"],rparm[dimnames(reg.proto)[[2]]])
                 },
                     rparm,
                     if (method=="F2") c(0,rparm[modes.used]) else NULL,
                     stop("unamed 'rparm' has wrong length")
                       )
      }
  }
  if (rparm.named&&is.na(rparm["intercept"]))
      rparm["intercept"] <- 0.0 # provide unstated intercept
  if (any(is.na(rparm)))
      stop(paste("rparm undefined for term:",reg.vec[is.na(rparm)][1]))
  if(nloc > 1) {
      need.prod <- as.numeric(ana.obj$loc.right[subset, z.used[ - nloc], drop = FALSE] > 
                              rep(z.used[-1], rep(n, nloc - 1)))
      lambda <- rep(0, nloc - 1) 
      for(i in 2:nloc)
          lambda[i - 1] <-
              prod(marker.distances[z.used[i - 1]:(z.used[i] - 1)])
  }
  else need.prod <- lambda <- 0
  if (is.null(ana.obj$casewt)) {
      casewt <- -1.0 #flag for no weights
  }
  else
  {
      if (length(ana.obj$casewt)!=length(subset))
          stop("length(casewt) must  == 1 or nrow(ana.obj$data)")
      else
          casewt <- ana.obj$casewt[subset]
  }
  if (setup.only) {
      res <-
          list("lapadj",
               csrtype=as.integer(crstype),
               nparm = as.integer(nparm),
               y = as.double(y),
               x = as.double(reg.proto),
               z = as.double(if (ncz>0) covar.matrix[,-1] else 0),
               rpar = as.double(rparm),
               tab =
               if (nrx>1) {
                   as.double(ana.obj$state.matrix[subset, z.used,  ])
               }
               else {
                   as.double(rep(1,n))
               } ,
               needProd = as.integer(need.prod),
               lambda = as.double(lambda),
               coefs = coefs,
               casewt=as.double(casewt),
               ln.sigma = ln.sigma,
               hkExact = double(1),
               hkApprox = double(1),
               llk = double(1),
               hess = double((np)^2),
               postApprox = double(1),
               iter = as.integer(iter),
               tol = as.double(tol),
               nem = as.integer(nem) )
      attr(res,"subset") <- subset
      attr(res,"reg.names") <- dimnames(reg.proto)[[2]]
      attr(res,"N")<-c(N=n.data,N.omit=n.data-n,N.used=n)
      res
  }
  else
  {
      z <- .C("lapadj",
              csrtype=as.integer(crstype),
              nparm = as.integer(nparm),
              y = as.double(y),
              x = as.double(reg.proto),
              z = as.double(if (ncz>0) covar.matrix[,-1] else 0),
              rpar = as.double(rparm),
              tab =
              if (nrx>1) {
                  as.double(ana.obj$state.matrix[subset, z.used,  ])
              }
              else {
                  as.double(rep(1,n))
              } ,
              needProd = as.integer(need.prod),
              lambda = as.double(lambda),
              coefs = coefs,
              casewt=as.double(casewt),
              ln.sigma = ln.sigma,
              hkExact = double(1),
              hkApprox = double(1),
              llk = double(1),
              hess = double((np)^2),
              postApprox = double(1),
              iter = as.integer(iter),
              tol = as.double(tol),
              nem = as.integer(nem),
              PACKAGE="bqtl")
      if ( (z$postApprox==0) && (z$hkApprox==0))
        adj <- 1.0
      else
        adj <- z$postApprox/z$hkApprox
      list(adj = adj, logpost = z$llk, parm =
           c(z$coefs[1: (nreg)], z$ln.sigma), posterior = z$postApprox,
           hk.approx = z$ hkApprox, hk.exact = z$hkExact, reg.vec = reg.vec,
           rparm = rparm, hess = if(return.hess) matrix(z$hess, nc = np) else NULL,
           iter = z$iter,N=c(N=n.data,N.omit=n.data-n,N.used=n),
           call=this.call)
  }
} 
