"map.dx"<-
    function(lambda = (1 - 2 * theta), theta = NULL, min.lambda = 0)
{
    if(length(lambda) && any(lambda < min.lambda)) {
        lambda <- pmax(lambda, min.lambda)
        warning("small value(s) for lambda corrected")
    }
    - log(lambda)/2
}
