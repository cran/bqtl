configs<-
    function(x,...,scope,method=NULL)
{
###
### configs will return a single variable name when x has length one
### and there are no additional arguments,
###  i.e. configs(x[1],scope=names(vars)) == list(names(vars)[x[1]])
###  additional args are concat'ed w "+" separators
### when dim(x)==NULL or dim(x)[[1]]==1 a vector is returned
###  i.e. configs(x,scope=names(vars)) == names(vars)[x]
### finally when dim(x)[[1]]>1 or dim(x$configs)[[1]]>1
### a vector is returned with elements paste(scope[x[,i]],collapse="+"
###
###  I parenthesize everything. S+3.4 needs this workaround 
###    ~(a):(u+v):(w) parses correctly, but not  ~a:(u+v):w

    if (missing(scope)) stop("missing scope arg - possible bqtl syntax error")
    if ( is.atomic(x) && (length(x) < 2) ) {
        res <- paste(scope[c(x,...)],collapse="+")
        return(paste("(",res,")",sep=""))
    }
    else {
        if (length(list(...))!=0)
            stop("only one arg allowed with vector or matrix")
    }
    
    if (is.element("configs",names(x)))
    {
        x <- x$configs
    }
    if (length( dm <- dim(x))>0)
    {
        dim(x) <- c(dm[1],prod(dm[-1]))
    }
    else
    {
        dim(x) <- c(1,length(x))
    }
    res <-    apply(x,2,function(x,dat) paste(dat[x],collapse="+"),dat=scope)
    paste("(",res,")",sep="")
}

