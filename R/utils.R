## mygene utility functions
library(Hmisc)
library(plyr)

.collapse <- function(...) {
    paste(unlist(list(...)), sep=",", collapse=",")
}

.transpose.nested.list <- function(li) {
    ## Assumes that inner names of each element are the same
    #if (length(li) == 0)
    #  return(li)
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

.pop <- function(list, item, default_value=NULL){
    if (is.null(list[[item]])){
        return(default_value)
    }
    else{
        value <- list[[item]]
        return(value)}
}

.unnest <- function(list) {
    while(any(vapply(list, is.list, TRUE))){
    list<-lapply(list, unlist, recursive=FALSE)
    return(list)
    }
}


.unnest.df <- function(df, recursive=TRUE) {
    reslist <-lapply(colnames(df), function(i) {
        if (is(df[[i]], "data.frame")) {
          if (recursive){
            df[[i]]<-.unnest.df(df[[i]], recursive=TRUE)
          }
            setNames(df[[i]], paste(i, colnames(df[[i]]), sep="."))
        } 
        else {
            df[i]
        }
    })
    res <- do.call(cbind, reslist)
    row.names(res) <- row.names(df)
    res
}

.json2df <- function(x){
    li <- lapply(x, fromJSON, flatten=TRUE)
    df <- rbind.fill(li)
    df
}

#before writing to TSV/CSV/xlsx
# .convert2csv<-function(df){
#     needpc <-sapply(df, is, "CharacterList")
#     df[needpc]<-lapply(df[needpc],rtracklayer:::pasteCollapse)
# }

.json.batch.collapse <- function(x){
    #stopifnot(all(grepl("^\\s*\\[.*\\]\\s*$", x, perl=TRUE)))
    x <- gsub(pattern="^\\s*\\[|\\]\\s*$", replacement="", x, perl=TRUE)
    x <- paste(x, collapse=",")
    paste("[", x, "]")
}

.uncollapse <- function(x, sep=",") {
    x <- as.character(unlist(x))
    unlist(strsplit(x, sep, fixed=TRUE))
}
