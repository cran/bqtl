"map.names" <-
    function(x,...)
{
    UseMethod("map.names")
}
"map.names.default"<-
    function(x,...)
{
    names(x,...)
}
"map.names.map.frame" <-
    function(x,...)
{
    c(x$marker.name[map.index(x,...)])
}
"map.names.analysis.object"<-
    function(x,...)
    map.names(x$map.frame,...)
"map.names.bqtl"<-
    function(x,...,ana.obj=NULL)
{
    if (exists("is.R") && is.R())
        tokens <- unlist(strsplit(x$reg.vec,":"))
    else
        tokens <-  unlist(lapply(x$reg.vec,strsplit,":"))
    covar.used <- tokens[grep("covar\\(",tokens)]
    
    if (is.null(ana.obj))
        ana.obj <- eval(match.call(bqtl,x$call)$ana.obj)
    
    
    reg.names.indx <- ana.obj$reg.names%in%tokens
    if (length(covar.used) != 0)
        reg.names.indx <-
            reg.names.indx | ( paste("covar(",ana.obj$reg.names,")",sep="")%in%
                              covar.used)
    res <-
        if (ana.obj$method == "F2")
            unique(
               {x<-dimnames(ana.obj$reg.names)[[2]];
                rep(x,rep(2,length(x)))[reg.names.indx]})
        else
            unique(dimnames(ana.obj$reg.names)[[2]][reg.names.indx])
    res
}
"map.names.bqtl.list"<-
    function(x,...)
{
    res <- list()
    for (i in seq(along=x))
        res[[i]]<-map.names.bqtl(x[[i]],...)
    res
}
"map.location"<-
    function(x,...)
    UseMethod("map.location")
"map.loc"<-map.location
"map.location.default"<-
    function(x,y,chromo=NULL,cM=NULL,map.names=NULL,...)
{
    if (missing(y)){
        sw.arg <-    1+2*is.null(chromo)+is.null(map.names)
        y <- switch(sw.arg,
                    ##both non-null
                    stop("cannot use both chromo= and map.name ="),
                    ## chromo
                    if (all(chromo%in%x$chr.num)){
                        if (is.null(cM))
                            x$chr.num %in% chromo
                        else
                            map.index(x,chromo,cM)
                    }
                    else
                    stop("bad chromo number specified"),
                    ## map.name
                    if (all(map.names%in%x$marker.name)){
                        x$marker.name%in%map.names}
                    else
                    stop(paste("unrecognized map.name:",
                               map.names[!map.names%in%x$marker.name])),
                    ## both null use all names
                    seq(nrow(x))
                    )
    }
    res <- if (mode(y)=="character")
        map.location(x,map.names=y)
    else
        x[y,c("chr.num","cM","marker.name")]
    class(res) <- c("map.location",class(x))
    res
}
"map.location.analysis.object"<-
    function(x,...)
    map.location(x$map.frame,...)
"map.location.bqtl"<-
    function(x,...,ana.obj=NULL)
{
    if (is.null(ana.obj))
        ana.object <- eval(match.call(bqtl,x$call)$ana.obj)
    mnames <-    map.names(x,ana.obj=ana.object)
    map.location(ana.object,mnames)
    
}
"map.location.bqtl.list"<-
    function(x,...)
{
    res <- list()
    for ( i in seq(along=x))
        res[[i]] <- map.location.bqtl(x[[i]],...)
    res
}
