"dom"<-
    function(x,...,scope=scope,method=method)
{
    if (missing(scope)) stop("missing scope arg - possible bqtl syntax error")
    new.scope <- scope[grep("^dom",scope)]
    if (length(new.scope)==0) stop("no 'dom' terms found")
    configs(x,...,scope=new.scope,method=method)
}
