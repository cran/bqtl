"bqtl"<-
    function(reg.formula,ana.obj=analysis.object, 
             scope = ana.obj$reg.names, expand.specials= NULL, ...)
{
    local.covar <-
        function (x, ..., scope = scope, method = method,
                  prefix = NULL, bq.spec=bqtl.specials) 
        {
            if (missing(scope)) 
                stop("missing scope arg - possible bqtl syntax error")
            x.arg <- match.call()$x
            if  ( length(x.arg) ==1 && deparse(x.arg) == "." )  # allow . shorthand 
                x <- seq(along=scope)
            if (is.call(x.arg)){
                if (deparse(x.arg[[1]])%in%bq.spec){
                    if (is.null(x.arg$scope)) x.arg$scope <- as.name("scope")
                    if (is.null(x.arg$method)) x.arg$method <- method
                    if (!missing(...)) stop("cannot use ... in this context")
                    paste("covar(",eval(x.arg),")",sep="")
                }
                else{ ## x.arg is c(1,2) or 7:8 or whatever
                    x <- eval(x.arg)
                    new.scope <- paste("covar(", scope, ")", sep = "")
                    locus(x, ..., scope = new.scope, method = method)
                }
            }
            else {
                deparse.x <- deparse(x.arg)
                if (deparse.x%in%scope){
                    if (!missing(...))
                        stop("cannot use ... args with named locus")
                    else
                        paste("covar(",deparse.x,")",sep="")
                }
                else {
                    new.scope <- paste("covar(", scope, ")", sep = "")
                    locus(x, ..., scope = new.scope, method = method)
                }
            }
        }
    
    this.call <- match.call()
    if (any(is.na(pmatch(c("loc.right","map.frame","state.matrix","method"),names(ana.obj)))))
        stop("ana.obj doesn't have required components")

    using.R <- exists("is.R")&&is.R()

### get terms.object, extract 'configs' terms,
    
    bqtl.specials <- c("configs","locus","add","dom","covar","acovar","dcovar")
    reg.terms <- terms(reg.formula,specials=bqtl.specials)
    reg.labels <- labels(reg.terms)
    reg.specials <- attr(reg.terms,"specials")  
    if (length(unlist(reg.specials))==0) { # no specials - pass thru to lapadj
        res <- 
            lapadj(reg.formula,ana.obj,...)
        res$call <- attr(res,"call")<- this.call
        class(res)<- "bqtl"
        return(res)
    }

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
    reg.plain <- term.rownames[ - reg.specials ]
    names(reg.plain) <- reg.plain
    response  <- term.rownames[1]
    reg.plain[response] <- paste(response, "~")
    
### let the specials expand themselves
    if (using.R)
        formals(local.covar)$bq.spec <- bqtl.specials # bind bqtl.specials explicitly
    else
        local.covar$bq.spec <- bqtl.specials
    pt.vars <- reg.specials + if (using.R) 1 else 0
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

    if (using.R ){
        term.list <- c(rspec,as.list(c(reg.plain,plus="+",colon=":")))
    }
    else { #argh! Splus 3.4
        rspec <- lapply(rspec,as.character)
        term.list <- c(rspec,as.list(c(reg.plain,plus="+",colon=":")))
        names(term.list) <- c(names(rspec),names(reg.plain),"plus","colon")
    }
    
### order is <var,conj,var,conj,...,conj,var>
    spec.col.order <-
        if (length(terms.conjuncs)==0)
            term.names
        else
            c(term.names,terms.conjuncs)[c(rep(c(0,n.join+1),n.join)+
                                           rep(1:n.join,rep(2,n.join)),n.join+1)]
####prepend response
    spec.col.order <- c(response, spec.col.order) 
    text.formula <- do.call("paste",term.list[spec.col.order])

### now process it all
    res <- list()
    for (i in seq(along=text.formula)) {
        this.formula <- eval(parse(text=text.formula[i]))
        res[[i]] <-
            lapadj(this.formula, ana.obj,...)
        class(res[[i]])<-"bqtl"
        res[[i]]$call <- substitute(bqtl(this.formula,ana.obj,...))
    }
    if (length(res)==1) {
        res <- res[[1]]
    }
    else
        class(res) <- "bqtl.list"

    attr(res,"call") <- substitute(this.call)
    res
}
