"twohk" <-
    function(varcov, ana.obj, ...)
{
    this.call <- match.call()
    if (!inherits(ana.obj, "analysis.object"))
        stop("second arg must be analysis.object")
    method <- ana.obj$method
    res <-
        if (method == "F2") {
            twohkf2(varcov, ana.obj, ...)
        }
        else {
            twohkbc1(varcov, ana.obj, ...)
        }
    res$call <- this.call
    res
}
