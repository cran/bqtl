"rhs.bqtl"<-
    function(reg.terms,ana.obj,bqtl.specials,local.covar,
             scope, expand.specials=NULL,method,...)
{
    reg.labels <- labels(reg.terms)
    reg.specials <- attr(reg.terms,"specials")  
    if (missing(method))
        method <- ana.obj$method
### work on attr(,"factors")
    reg.factor <- attr(reg.terms,"factors")
    term.rownames <- dimnames(reg.factor)[[1]]
    term.names <- term.rownames[row(reg.factor)[reg.factor>0]]
    terms.conjuncs <-
        if ( sum(reg.factor>0) <2 )
            NULL
        else
            ifelse( diff(col(reg.factor)[reg.factor>0])==0,"colon", "plus")
    n.join <- length(terms.conjuncs)
### flag rows with specials and those without
    reg.specials <- unlist(reg.specials[ !sapply(reg.specials,is.null) ])
    reg.plain <-
        if (is.null(reg.specials))
            term.rownames
        else
            term.rownames[ - reg.specials ]
    names(reg.plain) <- reg.plain
    
    
### let the specials expand themselves
    formals(local.covar)$bq.spec <- bqtl.specials # bind bqtl.specials explicitly
    
    pt.vars <- reg.specials +  1 
    if (length(reg.specials) != 0) {
        rspec <-
            lapply(attr(reg.terms,"variables")[pt.vars],
                   function(x,scope,method,covar) {
                       if ( !is.element("scope",names(x)) ) #typically use default
                           x$scope <- as.name("scope")
                       if ( !is.element("method",names(x)) ) #typically use default
                           x$method <- method
                       eval(x)
                   },
                   scope=scope,method=method,
                   covar=local.covar)
        names(rspec) <- term.rownames[ reg.specials ]
### used a common 'scope' and 'method' for all specials
        rspec.check <-
            sapply(rspec, function(x) any(x %in% c("(NA)","()")))
        if (any(rspec.check)){
            bad.terms <-
                paste(c("invalid term(s) in formula:",names(rspec)[rspec.check]),
                      collapse=" ")
            stop(bad.terms)
        }
        
### use all combinations of the expanded variables ?
        if (is.null(expand.specials))
            expand.specials <- 
                length(rspec)>1 && any(diff(range(sapply(rspec,length)))!=0)

        if (expand.specials)
            rspec <- do.call("expand.grid",rspec)
        else
            rspec <- lapply(rspec,paste,collapse="+")
    }
    else { # no specials
        rspec <- NULL
    }
    
    term.list <- c(rspec,as.list(c(reg.plain,plus="+",colon=":")))
    
### order is <var,conj,var,conj,...,conj,var>
    spec.col.order <-
        if (length(terms.conjuncs)==0)
            term.names
        else
            c(term.names,terms.conjuncs)[c(rep(c(0,n.join+1),n.join)+
                                           rep(1:n.join,rep(2,n.join)),n.join+1)]
    do.call("paste",term.list[spec.col.order])
}
