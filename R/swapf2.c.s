"swapf2" <-
    function(varcov,invars,rparm,nreps,ana.obj,
             locs=NULL,
             locs.prior = NULL,
             combo.prior = rep(1, 3)/3,tol=1e-10)
{
    this.call <- match.call()
    if (any(duplicated(invars)))
        stop("duplicate values of invars not allowed")
    nin<- length(invars)
    if (!missing(ana.obj)) {
        if (class(ana.obj) != "analysis.object")
            stop("arg \'ana.obj\' must be analysis object")
        if (missing(locs.prior)) {
            match.names <-
                match(dimnames(varcov$var.x)[[1]],ana.obj$reg.names,0)
            match.cols <- unique(col(ana.obj$reg.names)[match.names])
            if (any(match.names==0))
                stop("variable names in varcov not found in ana.obj")
            locs.prior <-
                ana.obj$map.frame[match.cols,"prior"]
            locs.prior <- locs.prior/sum(locs.prior)
        }
    }
    else {
        if (missing(locs.prior))
            locs.prior <- rep(1, ncol(varcov$var.x)/2)
    }

    if (missing(locs))
        locs <- rep(seq(ncol(varcov$var.x)/2), rep(2, ncol(varcov$var.x)/2))
    
    in.locs <- locs[invars]
    uniq.locs <- unique(in.locs)
    loc.tab <- table(in.locs)
    mode.offset <- ifelse(loc.tab[as.character(in.locs)]==2,2,(invars-1)%%2)
    mode.offset <- mode.offset[ match( uniq.locs, in.locs ) ]
    nlocs <- length(uniq.locs)
    nstep <- nlocs
    len.locs <- length(unique(c( locs) ))
    
    nopt <- len.locs*3
    optmax <- 2
    nrx <- nrow(varcov$var.x)
    if (length(combo.prior)!=3) stop("combo.prior must = 3")
    if (nrx!=length(rparm)) {
        if (length(rparm)<3)
            rparm <- rep(rparm,nrx/length(rparm))
        else
            stop(paste("rparm must have length 1, 2, or ",nrx))
    }
    if (length(locs.prior)!=len.locs)
        stop("length(locs.prior)!=length(unique(locs))")
    if (length(nreps)>1) stop("nreps should be an integer")
    optpri <- c(locs.prior)%o%c(combo.prior)
    optvars <- matrix(c(rbind(unique(c(locs))*2-2,-1),rbind(unique(c(locs))*2-1,-1),
                        locs*2-rep(2:1,nopt/3)),nrow=2)

    optcur <- uniq.locs + mode.offset*len.locs
    noptuse <- length(optcur)
    optblock <-
        rbind(rep(1:len.locs,3)+rep(c(len.locs,0,0),rep(len.locs,3)),
              rep(1:len.locs,3)+
              rep(c(2*len.locs,2*len.locs,len.locs),rep(len.locs,3)))-1
    useopt <- rep(1,nopt)
    useopt[c(optcur,optblock[,optcur]+1)] <- 0
    inpri <- 1
    xmax <- 2*(nlocs-1)
    zmax <- 2
    nmax <- xmax+zmax
    invars <- c(invars-1,rep(-1,nmax-nin))
    
    z <- .C("swapf2",
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
            alt.coef=double(nrx))

    res <- z[c("configs","posteriors","coefs","cond","marg","alt.marginal","alt.coef")]
    res$configs <- array(res$configs+1,c(nmax,nlocs,nreps))
    dim(res$coefs) <- c(nmax,nlocs,nreps)
    dim(res$marg) <- dim(res$cond) <- dim(res$posteriors) <- c(nlocs,nreps)
    res$nloc <- len.locs
    res[["call"]] <- this.call
    class(res) <- "swap"
    res
}
