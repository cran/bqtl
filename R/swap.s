"swap" <-
    function(varcov, invars, rparm, nreps, ana.obj, ...)
{
    this.call <- match.call()
    if (class(ana.obj) != "analysis.object")
        stop("arg 5 must be analysis.object")
    method <- ana.obj$method
    res <-
        if (method == "F2") {
            swapf2(varcov, invars, rparm, nreps, ana.obj, ...)
        }
        else {
            swapbc1(varcov, invars, rparm, nreps, ana.obj, ...)
        }
    res$call <- this.call
    res
}
