if (!exists("is.element")) {
    "is.element" <- function (el, set) match(el, set, 0) > 0
    "%in%"<-is.element
}
