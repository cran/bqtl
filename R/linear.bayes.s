"linear.bayes" <-
    function(x, ana.obj, partial=NULL,rparm,specs,
             scope,subset,casewt,...){
        x.call <- match.call(expand.dots=TRUE)
        if (class(x)!="varcov"){
            vc.call <- x.call
            vc.call[[1]] <- as.name("varcov")
            OK.args <-
                is.element(names(vc.call),
                           c("subset","casewt","ana.obj","partial",
                             "scope","x"))
            OK.args[1] <- TRUE
            vc.call[!OK.args] <- NULL
            x <- eval(vc.call)
        }
        else {
            vc.call <- x$call
        }
        default.specs <-
            list(
                 gene.number=c(1,2,3,4,5),
                 burn.in=1,
                 n.cycles=c(0,200,200,100,100))
        
        if (missing(specs)){              
            specs <- default.specs
        }
        else {
            named <- pmatch( names(specs), names(default.specs), 0)
            if (any(named == 0)) stop("illegal name in specs arg")
            names(specs) <- names(default.specs)[named]
            specs <- c(specs, default.specs[-named] )
        }
        do.hk <- any(specs$gene.number<3)
        if (missing(rparm)) {
            warning("default: rparm  <- 0  may be unstable")
            rparm <- 0
        }
        if (do.hk)
            hk <- twohk(x,ana.obj,rparm=rparm,...)
        else
            hk <- NULL
        nums.2.sample <- sort(specs$gene.number[specs$n.cycles!=0])
        contig.nums <- all(diff(nums.2.sample)==1)
        if (!contig.nums) stop("specs$gene.number must be contiguous")
        var.indx <-
            if (ana.obj$method=="F2")
                seq(from=1,by=2,length=ncol(x$var.x)/2)
            else
                seq(ncol(x$var.x))
        if (specs$burn.in > 0){
            starts <- vector("list",max(nums.2.sample))
            for (i in nums.2.sample) {
                tmp <-  swap(x,sample(var.indx,i),
                             rparm,1,ana.obj=ana.obj,...)$configs[,i,1]
                starts[[i]] <- tmp[tmp != 0]
            }
            
        }
        swaps <- vector("list",max(nums.2.sample))
        for (i in nums.2.sample){
            j <- match(i,specs$gene.number)
            cycles <- specs$n.cycles[[j]]
            if (cycles != 0)
                swaps[[i]] <-
                    swap(x,starts[[i]], rparm,cycles, ana.obj=ana.obj,...)
        }
        swaps.OK <- sapply(swaps,function(x) 
                           if (length(unlist(x))==0)
                           TRUE          # NULL is OK
                           else {
                               call.indx <- match("call",names(x),0)
                               !any(is.na(unlist( x[-call.indx] )))})
        smry <- vector("list",length(swaps))
        
        for (i in seq(along=swaps))
            if (swaps.OK[[i]]) smry[[i]] <- summary(swaps[[i]])

        if (any(!swaps.OK)){
            warning("attempt to summarize samples failed. rparm = 0 ??")
            odds <- loc.posterior <- coefs <- NULL          
        }
        else{
            swap.odds <- sapply(smry,function(x) 1/x$ratio$mean)
            if (do.hk){
                odds <-
                    c(sum(hk$loc.1),
                      cumprod(c(sum(hk$loc.2),
                                swap.odds[nums.2.sample[nums.2.sample>2]])))
                post.pr <- odds/sum(odds)
                if (ana.obj$method=="F2") {
                    loc.posterior <-
                        c(hk$loc.1%*%rep( 1/sum(hk$loc.1), 3)) * post.pr[1] +
                            2*c(hk$loc.2%*%rep( 1/sum(hk$loc.2), 3)) *
                                post.pr[2] +
                                    sapply(smry[-(1:2)],
                                           "[[", "loc.posterior" ) %*%
                                               post.pr[-(1:2)]
                }
                else {
                    loc.posterior <-
                        hk$loc.1 / sum(hk$loc.1) * post.pr[1] +
                            2 * hk$loc.2 / sum(hk$loc.2) * post.pr[2] +
                                sapply(smry[-(1:2)],
                                       "[[", "loc.posterior" ) %*%
                                           post.pr[-(1:2)]
                }
                coefs <-
                    cbind(c(hk$coefs.1),c(hk$coefs.2),
                          sapply(smry[-(1:2)],"[[","coefs")) %*% post.pr
            }
            else {
                swap.nums <- nums.2.sample[nums.2.sample>2]
                odds <-
                    cumprod(c(1,swap.odds[ swap.nums[-1] ]))
                post.pr <- odds/sum(odds)
                loc.posterior <-
                    sapply(smry[swap.nums],"[[", "loc.posterior" ) %*% post.pr
                coefs <- sapply(smry[swap.nums],"[[","coefs") %*% post.pr
            }
            
        }
        res <- list(hk=hk, swaps=swaps, smry=smry, odds=odds,
                    coefficients=coefs, loc.posterior=loc.posterior,
                    specs=specs,varcov.call=vc.call,call=x.call)
        class(res) <- "linear.bayes"
        res
    }

