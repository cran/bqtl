"twohkbc1"<-
    function(varcov, ana.obj, rparm=0, locs = NULL ,
             locs.prior = NULL )
{
### step thru all locations and compute the marginal and posterior
### probability for each location
    this.call <- sys.call()
    if (!missing(ana.obj)) {
        if (class(ana.obj) != "analysis.object")
            stop("arg \'ana.obj\' must be analysis object")
        if (missing(locs.prior)) {
            match.names <-
                match(dimnames(varcov$var.x)[[1]],ana.obj$reg.names,0)
            if (any(match.names==0))
                stop("variable names in varcov not found in ana.obj")
            locs.prior <-
                ana.obj$map.frame[match.names,"prior"]
        }
    }
    else {
        if (missing(locs.prior))
            locs.prior <- rep(1, ncol(varcov$var.x))
    }
    if (missing(locs)) 
        locs <- seq(ncol(varcov$var.x))
    nlocs <- length(locs)
    nin <- 1
    nopt <-  nlocs
    optmax <- 1
    nmax <- 2
    nrx <- nrow(varcov$var.x)
    optpri <- c(locs.prior)
    optvars <- locs - 1
    marg <- post <- cond <- rep(0, nlocs)
    coefs <- rep(0, nlocs)
    useopt <- rep(1,nlocs)
    if (nrx!=length(rparm)) {
        if (length(rparm)==1)
            rparm <- rep(rparm,nrx)
        else
            stop(paste("rparm must have length 1 or ",nrx))
    }
    if (nlocs != length(locs.prior))
        stop("lengths of locs.prior and locs not equal")
    
    res <- .C("twohkbc1",
              varx = as.double(varcov$var.x),
              covxy = as.double(varcov$cov.xy),
              vary = as.double(varcov$var.y),
              df = as.double(varcov$df),
              rparm = as.double(rparm),
              optpri = as.double(optpri),
              nrx=as.integer(nrx),
              optvars=as.integer(optvars) ,
              nopt = as.integer(nopt),
              optmax = as.integer(optmax),
              useopt=as.integer(useopt),
              pvt = as.integer(1:nmax),
              rank = integer(1),
              wrksp = double(2*nmax),
              gama = double(nin*optmax),
              bee = double(nin),
              xx = double(nin*nin),
              xy = double(nin),
              zz = double(optmax*optmax),
              zy = double(optmax),
              xz = double(nin*optmax),
              beta =double(optmax),
              posterior = double(nopt),
              loc.2 =as.double(marg),
              loc.1 = as.double(cond),
              coefwk = double(nmax*nopt),
              coefs.2 = as.double(coefs),
              coefs.1=as.double(coefs),
              qraux =double(nmax),
              zraux =double(nmax),
              zrank=integer(1),
              tol=as.double(1e-10))[c("loc.2","loc.1","coefs.2","coefs.1")]
    dim(res$loc.2) <- c(nlocs,1)
    dim(res$loc.1) <- c(nlocs,1)
    dim(res$coefs.2) <- c(1,nlocs)
    dim(res$coefs.1) <- c(1,nlocs)
    res[["call"]] <- this.call
    res
}
