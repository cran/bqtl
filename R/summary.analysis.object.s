"summary.analysis.object" <-
    function(x)
{
    mf <- summary(x$map.frame)
    norg <- nrow(x$data)
    method <- x$method
    origin<-x$call
    pheno.cov <- setdiff(names(x$data),x$reg.names)
    list(map=mf,N=norg,method=method,variable.names=pheno.cov,
         mode.mat=x$mode.mat,call=origin)
}

