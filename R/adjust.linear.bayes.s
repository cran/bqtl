"adjust.linear.bayes" <-
    function(lbo,ana.obj=lbo$call$ana.obj,...)
{
    this.call <- match.call(expand.dots=TRUE)
    specs <- lbo$specs

### get twohk and swap args from lbo
    cl.1 <- lbo$hk$call
    cl.1[[1]]<-cl.1$varcov<-cl.1$ana.obj<-NULL
    cl.2 <- lbo$varcov.call
    cl.2[[1]]<-cl.2$x<-cl.2$ana.obj<-NULL

### this builds up the call to bqtl using args from
### the original call to linear.bayes. It figures out which to bother
### with by inspecting the calls to varcov and twohk as invoked by
### linear.bayes      
    one.gene.expr <- call("bqtl",reg.formula=lbo$varcov.call$x,
                          ana.obj=ana.obj)
### are there extra args to be added ??
    extra.args <- unique(c(names(cl.1),names(cl.2)))
    if (length(extra.args)!=0){
        dummy.call <- expression("$<-"(one.gene.expr,arg.2,arg.3))[[1]]
        for ( i in extra.args ){
            dummy.call[[3]] <- dummy.call[[4]] <- i
            eval(dummy.call)
            one.gene.expr[[i]]<-lbo$call[[i]]
        }
    }
    
### build call to summary.adj
    one.gene.smry.expr <-
        call("summary.adj",adj.obj=as.name("one.gene.models"),
             n.loc=1,coef.znames="",mode.names=NULL,imp.denom="")
    one.gene.smry.expr$coef.znames <- call("$",ana.obj,"reg.names")
    mf.expr <- call("$",ana.obj,"map.frame")
    prior.expr <- call("$",mf.expr,"prior")
    one.gene.smry.expr$imp.denom <- call("/",1,prior.expr)
    
    if ( any(specs$gene.number == 1) ) {
        one.gene.models <- eval(one.gene.expr)
        one.gene.smry <- eval(one.gene.smry.expr)
    }
    else {
        one.gene.smry <- NULL
    }
### now do multi-gene models
    swap.gene.expr <- one.gene.expr
    swaps.expr <-
        expression( configs( unique.config(lbo$swaps[[1]])$uniq ) )[[1]]
    swaps.expr[[2]][[2]][[2]][[2]][[2]] <- as.name(this.call$lbo)
    n.gene.smry.expr <- one.gene.smry.expr
    n.gene.smry.expr$adj.obj <- as.name("swaps.adj")
    n.gene.smry.expr$imp.denom <- NULL
    n.gene <- vector("list",length(lbo$swaps))
    n.swap <- specs$gene.number[specs$n.cycles != 0]

    for( i in n.swap ){
        ## y ~ configs( unique.config(lbo$swaps[[  i  ]])$uniq )
        swaps.expr[[2]][[2]][[2]][[3]] <- i
        swap.gene.expr$reg.formula[[3]] <- swaps.expr
        swaps.adj <- eval(swap.gene.expr)  # do it
        ## now summarize
        n.gene.smry.expr$n.loc <- i
        n.gene.smry.expr$swap.obj <- swaps.expr[[2]][[2]][[2]]
        n.gene[[i]] <- c( eval(n.gene.smry.expr), list(swaps=swaps.adj))
    }
    if ( is.null(lbo$odds) ){
        odds <- loc.posterior <- coefs <- NULL
    }
    else
        if ( any(specs$gene.number == 1) ){
            odds <-
                lbo$odds*cumprod(c(one.gene.smry$adj,
                                   sapply(n.gene[n.swap],"[[","adj")))
            post.pr <- odds/sum(odds)
            loc.posterior <-
                cbind(one.gene.smry$loc,sapply(n.gene[n.swap],"[[","loc")) %*%
                    ( post.pr* c(1,n.swap))
            coefs <-
                cbind(one.gene.smry$coef,sapply(n.gene[n.swap],"[[","coef")) %*%
                    post.pr
        }
        else {
            odds <-
                lbo$odds*cumprod( sapply(n.gene[n.swap],"[[","adj") )/
                    n.gene[n.swap][[1]]$adj  #take first as baseline for odds
            post.pr <- odds/sum(odds)
            loc.posterior <-
                sapply(n.gene[n.swap],"[[","loc") %*%
                    ( post.pr * n.swap )
            coefs <-
                sapply(n.gene[n.swap],"[[","coef") %*%
                    post.pr
        } 
    
    res <- list(odds=odds,loc.posterior=loc.posterior,coefficients=coefs,
                one.gene.adj=one.gene.smry,n.gene.adj=n.gene,
                call=this.call)
    class(res) <- "adjust.linear.bayes"
    
    res
}
