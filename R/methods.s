"logpost"<-
    function(x,...)
    UseMethod("logpost")
"loglik"<-
    function(x,...)
    UseMethod("loglik")


"logpost.bqtl" <-
    function(bqtl.obj,...)
{
    bqtl.obj$logpost
}

"logpost.default"<-logpost.bqtl
"loglik.bqtl"<-logpost.bqtl
"loglik.default"<-logpost.bqtl

"logpost.bqtl.list"<-
    function(bqtl.list,...)
{
    res <- sapply(bqtl.list,logpost.bqtl)
    res
}

"loglik.bqtl.list"<-logpost.bqtl.list

"coef.bqtl"<-
    function(bqtl.obj,...)
{
    ncoefs <- length(bqtl.obj$parm)-1
    coefs <- bqtl.obj$parm[1:ncoefs]
    names(coefs) <- c("Intercept",bqtl.obj$reg.vec)
    coefs
}

### method dispatch for coef takes care of this: "coefficients.bqtl"<-coef.bqtl

"coef.bqtl.list"<-
    function(bqtl.list,...)
{
    res <- sapply(bqtl.list,coef.bqtl)
    res
}

### method dispatch for coef takes care of this: "coefficients.bqtl.list"<-coef.bqtl.list

"summary.bqtl"<-
    function(bqtl.obj,...)
{
    nparm <-    length(bqtl.obj$parm)
    coefs <-
        if (is.null(bqtl.obj$hess))
            coef(bqtl.obj)
        else
        {
            vals <- coef(bqtl.obj)
            hess.qr <- qr(-bqtl.obj$hess)
            if (hess.qr$rank < nparm){
                use.coef <-  vals != 0
                hrc <- c(use.coef,TRUE)
                std.err <- ifelse( use.coef, 0, NA)
                std.err[use.coef]<-
                    sqrt(diag(solve(-bqtl.obj$hess[hrc,hrc]))[-hess.qr$rank])
            }
            else{
                std.err <- sqrt(diag(solve(hess.qr))[-nparm])
            }
            t.stats <- vals/std.err
            p.vals <- pnorm(-abs(t.stats))*2
            cbind(Value=vals,Std.Err=std.err,t.statistic=t.stats,p.value=p.vals)
        }
    
    loglik <- loglik(bqtl.obj)
    sigma <- exp(bqtl.obj$parm[nparm])
    list(coefficients=coefs,loglik=loglik,std.res=sigma,N=bqtl.obj$N)
}

"posterior"<-
    function(x,...)
    UseMethod("posterior")

"posterior.bqtl" <-
    function(bqtl.obj,...)
{
    bqtl.obj$posterior
}

"posterior.default" <-posterior.bqtl

"posterior.bqtl.list"<-
    function(bqtl.list)
{
    res <- sapply(bqtl.list,posterior.bqtl)
    res
}

