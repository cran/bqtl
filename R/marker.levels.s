### convenience functions for defining marker levels
###

"f2.levels" <-
    function(AA="AA", Aa="Aa", aa="aa", not.aa="A-", not.AA="a-",
             miss.val="--")
{
    this.call<-match.call()
    formal.names <- names(formals(eval(this.call[[1]])))
    res <- c(AA,Aa,aa,not.aa,not.AA,miss.val)
    names(res) <- formal.names
    res
}

"bc1.levels" <-
    function(AA="AA", Aa="Aa", miss.val="--")
{
    this.call<-match.call()
    formal.names <- names(formals(eval(this.call[[1]])))
    res <- c(AA, Aa, miss.val )
    names(res) <- formal.names
    c(res[1:2],rep("nil",3),res[3])
}

"ri.levels" <-
    function(AA="AA", aa="aa", miss.val="--")
{
    this.call<-match.call()
    formal.names <- names(formals(eval(this.call[[1]])))
    res <- c(AA, aa, miss.val )
    names(res) <- formal.names
    c(res[1:2],rep("nil",3),res[3])
}

