locus<-
    function(x,...,scope=scope,method=NULL,chromo=NULL,cM=NULL,ana.obj=NULL)
{
###
###   
### When method != "F2", locus will return a single variable name when x
### has length one, there are no additional arguments. When method ==
### "F2", then there will be two values returned if x has length one,
### hopefully the additive and dominance terms of that locus.
### 
### additional args are concat-ed w "+" separators when dim(x)==NULL or
### dim(x)[[1]]==1 a vector is returned.
### 
### I parenthesize everything. S+3.4 needs this workaround 
###  ~(a):(u+v):(w) parses correctly, but not  ~a:(u+v):w
###
  if (missing(scope)) stop("missing scope arg - possible bqtl syntax error")
  dots <- list(...)
  if (is.null(chromo)&& length(names(dots))!=0 &&
      any(cloc <- which(1==pmatch(names(dots),"chromo",0)))){
    chromo <- dots[[cloc]]
    if (length(dots)>1)
      stop("... not allowed")
    else
      dots<-NULL
  }
  if (!is.null(chromo)){
    if (!missing(x)) stop("using both x and chromo args not allowed")
    x <- if (is.null(cM))
        map.index(ana.obj,chromo=chromo)
    else
        map.index(ana.obj,chromo=chromo,cM=cM)
  }

  x.call <- match.call()$x
  if  ( length(x.call) ==1 && deparse(x.call) == "all" ){  # all loci ? 
    x <- seq(along=scope)
    if (method=="F2") dim(x) <- c(2,length(x)%/%2)
    if (length(dots) != 0) stop("... not allowed")
    return(configs(x,scope=scope))
  } #else
  
  if (length(x)>1 && length(dots)!=0)
    stop("only one arg allowed with vector or matrix")
  if (method == "F2") {
    if (length(x)==1) {
      x.1 <- 2*x-1
      y <- 2*x
      if ( length(dots)==0 )
        dots <- list(y)
      else
        dots <- c(list(y),lapply(dots,function(x) c(2*x-1,2*x)))
      configs.args <-  c(list(x=x.1),dots,scope=list(scope))
    }
        else {
            x.1 <- 2*x-1
            y <- 2*x
            if ( length( dm <- dim(x) )>0 && dm[1]>1 ) {
                x.1 <- aperm(array(c(x.1,y),c(dm[1],prod(dm[-1]),2)),c(3,1,2))
                dim(x.1) <- c(dm[1]*2,prod(dm[-1]))
            }
            else {
                x.1 <- rbind(x.1,y)
            }
            configs.args <-  c(list(x.1),scope=list(scope))
        }
    do.call("configs",configs.args)
  }
  else {
    configs.args <- c(list(x),dots,scope=list(scope))
    do.call("configs",configs.args)
  }
  
}

