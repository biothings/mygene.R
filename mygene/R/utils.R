library(Hmisc)

.collapse <- function(...) {
    paste(unlist(list(...)), sep=",", collapse=",")
}

.transpose.nested.list <- function(li) {
    ## Assumes that inner names of each element are the same
    inner.i <- seq_along(li[[1]])
    res <- lapply(inner.i, function(i) lapply(li, `[[`, i))
    names(res) <- names(li[[1]])
    res
}

.splitBySize <- function(x, maxsize) {
    n <- length(x)
    num.chunks <- ceiling(n / maxsize)
    f <- cut2(1:n, g=num.chunks)
    unname(split(x, f))
}

.df2DF <- function(df) {
    DF <- DataFrame(df, check.names=FALSE)
    isli <- sapply(df, is.list)
    DF[isli] <- lapply(df[isli], as, "List")
    DF
}
