"twohkf2"<-
    function(varcov, ana.obj, rparm=0, locs = NULL , locs.prior = NULL ,
             combo.prior = rep(1, 3)/3)
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
            match.cols <- unique(col(ana.obj$reg.names)[match.names])
            if (any(match.names==0))
                stop("variable names in varcov not found in ana.obj")
            locs.prior <-
                ana.obj$map.frame[match.cols,"prior"]
        }
    }
    else {
        if (missing(locs.prior))
            locs.prior <- rep(1, ncol(varcov$var.x)/2)
    }

    if (missing(locs))
        locs <- rep(seq(ncol(varcov$var.x)/2), rep(2, ncol(varcov$var.x)/2))
    
    uniq.locs <- unique(locs)
    nlocs <- length(uniq.locs)
    nin <- 2
    nopt <-  nlocs*3
    optmax <- 2
    nmax <- 4
    nrx <- nrow(varcov$var.x)
    if (length(combo.prior)!=3) stop("combo.prior must = 3")
    if (nrx!=length(rparm)) {
        if (length(rparm)<3)
            rparm <- rep(rparm,nrx/length(rparm))
        else
            stop(paste("rparm must have length 1, 2, or ",nrx))
    }
    if (length(locs.prior)!=nlocs)
        stop("length(locs.prior)!=length(unique(locs))")
    
    optpri <- c(locs.prior)%o%c(combo.prior)
    optvars <- matrix(c(rbind(seq(from=1,by=2,length=nlocs),0),
                        rbind(seq(from=2,by=2,length=nlocs),0),
                        seq(from=1,by=1,length=2*nlocs))-1,nrow=2)
    marg <- post <- cond <- rep(0, nlocs*3)
    coefs <- rep(0, 2*nlocs)
    useopt <- rep(1,nlocs*3)
    
    res <- .C("twohkf2",
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
    dim(res$loc.2) <- c(nlocs,3)
    dim(res$loc.1) <- c(nlocs,3)
    dim(res$coefs.2) <- c(2,nlocs)
    dim(res$coefs.1) <- c(2,nlocs)
    res[["call"]] <- this.call
    res
}
