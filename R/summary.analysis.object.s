"summary.analysis.object" <-
    function(object,...)
{
    mf <- summary(object$map.frame)
    norg <- nrow(object$data)
    method <- object$method
    origin<-object$call
    pheno.cov <- setdiff(names(object$data),object$reg.names)
    list(map=mf,N=norg,method=method,variable.names=pheno.cov,
         mode.mat=object$mode.mat,call=origin)
}

