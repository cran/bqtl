"swapbc1" <-
    function(varcov,invars,rparm,nreps,ana.obj,
             locs=NULL,
             locs.prior = NULL,
             tol=1e-10)
{
    this.call <- match.call()
    if (any(duplicated(invars)))
        stop("duplicate values of invars not allowed")
    nin <- length(invars)

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

            locs.prior <- locs.prior/sum(locs.prior)
        }
    }
    else {
        if (missing(locs.prior))
            locs.prior <- rep(1, ncol(varcov$var.x))
    }
    if (missing(locs)) 
        locs <- seq(ncol(varcov$var.x))

    in.locs <- unique(locs[invars])
    nlocs <- length(in.locs)
    nstep <- nlocs
    
    len.locs <- length(unique(c( locs ) ))
    
    nopt <- len.locs
    optmax <- 1
    nrx <- nrow(varcov$var.x)
    if (nrx!=length(rparm)) {
        if (length(rparm)==1)
            rparm <- rep(rparm,nrx)
        else
            stop(paste("rparm must have length 1 or ",nrx))
    }
    if (length(locs.prior)!=len.locs)
        stop("length(locs.prior)!=length(unique(locs))")
    if (length(nreps)>1) stop("nreps should be an integer")
    optpri <- c(locs.prior)
    optvars <- locs - 1
    optcur <- in.locs
    noptuse <- length(optcur)

    optblock <- (1:nopt)-1
    useopt <- rep(1,nopt)
    useopt[c(optcur,optblock[optcur]+1)] <- 0
    inpri <- 1
    xmax <- (nlocs-1)
    zmax <- 1
    nmax <- xmax+zmax
    invars <- c(invars-1)
    z <- .C("swapbc1",
            nreps=as.integer(nreps),
            nstep=as.integer(nstep),
            varx=as.double(varcov$var.x),
            covxy=as.double(varcov$cov.xy),
            vary=as.double(varcov$var.y),
            df=as.double(varcov$df),
            rparm=as.double(rparm),
            optpri=as.double(optpri),
            inpri=as.double(inpri),
            nrx=as.integer(nrx),
            noptuse=as.integer(noptuse),
            optblock=as.integer(optblock),
            optcur=as.integer(optcur-1),
            invars=as.integer(invars),
            nin=as.integer(nin),
            optvars=as.integer(optvars),
            nopt=as.integer(nopt),
            optmax=as.integer(optmax),
            useopt=as.integer(useopt),
            locs=as.integer(locs-1),
            pvt=as.integer(1:nmax),
            rank=integer(1),
            wrksp=double(2*nmax),
            gama=double(zmax*xmax),
            bee=double(xmax),
            xx=double(xmax^2),
            xy=double(xmax),
            zz=double(zmax^2),
            zy=double(zmax),
            xz=double(xmax*zmax),
            beta=double(xmax),
            posteriors=double(nreps*nstep),
            marg=double(nreps*nstep),
            cond=double(nreps*nstep),
            coefs=double(nreps*nstep*nmax),
            configs=integer(nreps*nstep*nmax),
            qraux=double(xmax),
            zraux=double(zmax),
            zrank=integer(1),
            tol=as.double(tol),
            postwk=double(nopt),
            coefwk=double(nopt*nmax),
            alt.marginal=double(nopt),
            alt.coef=double(nrx),
            PACKAGE="bqtl" )
    res <- z[c("configs","posteriors","coefs","cond","marg",
               "alt.marginal","alt.coef")]
    res$configs <- array(res$configs+1,c(nmax,nlocs,nreps))
    dim(res$coefs) <- c(nmax,nlocs,nreps)
    dim(res$marg) <- dim(res$cond) <- dim(res$posteriors) <- c(nlocs,nreps)
    res$nloc <- nopt
    res[["call"]] <- this.call
    class(res)<-"swap"
    res
}
