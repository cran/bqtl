"logpost"<-
    function(x,...)
    UseMethod("logpost")
"loglik"<-
    function(x,...)
    UseMethod("loglik")


"logpost.bqtl" <-
    function(x,...)
{
    x$logpost
}

"logpost.default"<-logpost.bqtl
"loglik.bqtl"<-logpost.bqtl
"loglik.default"<-logpost.bqtl

"logpost.bqtl.list"<-
    function(x,...)
{
    res <- sapply(x,logpost.bqtl)
    res
}

"loglik.bqtl.list"<-logpost.bqtl.list

"coef.bqtl"<-
    function(object,...)
{
    ncoefs <- length(object$parm)-1
    coefs <- object$parm[1:ncoefs]
    names(coefs) <- c("Intercept",object$reg.vec)
    coefs
}

### method dispatch for coef takes care of this: "coefficients.bqtl"<-coef.bqtl

"coef.bqtl.list"<-
    function(object,...)
{
    res <- sapply(object,coef.bqtl)
    res
}

### method dispatch for coef takes care of this: "coefficients.bqtl.list"<-coef.bqtl.list

"summary.bqtl"<-
    function(object,...)
{
    nparm <-    length(object$parm)
    coefs <-
        if (is.null(object$hess))
            coef(object)
        else
        {
            vals <- coef(object)
            hess.qr <- qr(-object$hess)
            if (hess.qr$rank < nparm){
                use.coef <-  vals != 0
                hrc <- c(use.coef,TRUE)
                std.err <- ifelse( use.coef, 0, NA)
                std.err[use.coef]<-
                    sqrt(diag(solve(-object$hess[hrc,hrc]))[-hess.qr$rank])
            }
            else{
                std.err <- sqrt(diag(solve(hess.qr))[-nparm])
            }
            t.stats <- vals/std.err
            p.vals <- pnorm(-abs(t.stats))*2
            cbind(Value=vals,Std.Err=std.err,t.statistic=t.stats,p.value=p.vals)
        }
    
    loglik <- loglik(object)
    sigma <- exp(object$parm[nparm])
    list(coefficients=coefs,loglik=loglik,std.res=sigma,N=object$N)
}

"posterior"<-
    function(x,...)
    UseMethod("posterior")

"posterior.bqtl" <-
    function(x,...)
{
    x$posterior
}

"posterior.default" <-posterior.bqtl

"posterior.bqtl.list"<-
    function(x,...)
{
    res <- sapply(x,posterior.bqtl)
    res
}

