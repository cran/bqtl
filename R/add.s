"add"<-
    function(x,...,scope,method=method)
{
    if (missing(scope)) stop("missing scope arg - possible bqtl syntax error")
    new.scope <- scope[grep("^add",scope)]
    if (length(new.scope)==0) stop("no 'add' terms found")
    configs(x,...,scope=new.scope,method=method)
}
