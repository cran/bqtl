"varcov"<-
    function(x,ana.obj,partial=NULL,scope=ana.obj$reg.names,...)
{
    this.call <- match.call(expand=TRUE)
    using.R <- exists("is.R")&&is.R()
    bqtl.specials <- c("configs","locus") # dont use: add or dom so
                                        # swap, hk wont jam up    
    reg.terms <- terms(x,specials=bqtl.specials)
    if (any(attr(reg.terms,"order")>1))
        stop("cannot use interactions in this formula")

    covs <-
        rhs.bqtl(reg.terms,ana.obj,bqtl.specials,NULL,scope,
                 expand.specials=FALSE)
    response <-
        dimnames(attr(reg.terms,"factors"))[[1]][attr(reg.terms,"response")]
    if (!is.null(partial)) {
        p.terms <- terms(partial,specials=bqtl.specials)
        if (any(attr(p.terms,"order")>1))
            stop("cannot use interactions in this formula")
        cntl <-
            rhs.bqtl(p.terms,ana.obj,bqtl.specials,NULL,scope,
                     expand.specials=FALSE)
        dummy.formula <- eval(parse(text=paste(response,"~",covs,"+",cntl)))
    }
    else {
        cntl <- NULL
        dummy.formula <- eval(parse(text=paste(response,"~",covs)))
    }
    dummy.terms <- terms(dummy.formula)
    mf <- model.frame(dummy.terms,ana.obj$data, na.action=na.omit)
    if (is.null(partial)) {
        x <- model.matrix(dummy.terms,mf)[,-1]
        y <- model.extract(mf,"response")
        make.varcov(x,y,...)
    }
    else
    {
        
        cntl.out <- eval(parse(text=paste("~ . -", cntl)))
        x.formula <- update(dummy.formula,
                            eval(parse(text=paste("~ . -", cntl)))  )
        x <- model.matrix(terms(x.formula),mf)[,-1]
        y <- model.extract(mf,"response")
        z <- model.matrix(terms(eval(parse(text=paste("~", cntl)))),mf)[,-1,drop=F]
        y[] <- lsfit(z,y)$resid
        x[] <- lsfit(z,x)$resid
        res <- make.varcov(x,y,...)
        r <- ncol(z)
        res$df <- res$df - r
        res$call <- this.call
        ## further adjustments unneeded since posterior is scaled by variance
        res
    }
    
}
