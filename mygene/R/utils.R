library(Hmisc)
library(xlsx)

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

.pop <- function(list, item, default_value=NULL){
    if (is.null(list[[item]])){
        return(default_value)}
    else{
        value <<- list[[item]]
        return(value)}
}

.unnest <- function(list) {
  while(any(vapply(list, is.list, T))){
  list<-lapply(list, unlist, recursive=FALSE)
  return(list)}
}

.unnest.df<-function(df){
    outdf <-jsonlite:::simplify(df)
    .df2DF(data.frame(as.list(outdf), check.names=FALSE))}

#before writing to TSV/CSV/xlsx
.convert2csv<-function(df){
    needpc <-sapply(df, is, "CharacterList")
    df[needpc]<-lapply(df[needpc],rtracklayer:::pasteCollapse)
}
#suggested
#write.xlsx(df, "out.xlsx", row.names=FALSE)

# uncollapse <- function(x, sep=",") {
#     unlist(strsplit(x, sep, fixed=TRUE))

# }







