library(IRanges)
library(httr)
library(jsonlite)

version <- '0.3'
MyGene <- setClass("MyGene",
    slots=list(base.url="character", delay="numeric", step="numeric", version="character", verbose="logical", debug="logical"),
    prototype=list(base.url="http://mygene.info/v2", delay=1, step=1000, version=version, verbose=TRUE, debug=FALSE))

validMyGeneObject <- function(object) {
    errors <- character(0)
    for (sn in c("base.url", "delay", "step")) {
        if (length(slot(object, sn)) != 1)
            errors <- c(errors, sprintf("Slot %s must have length 1", sn))
    }
    if (length(errors) > 0)
        errors
    else
        TRUE
}

setValidity("MyGene", validMyGeneObject)

.return.as<-function(gene_obj, return.as = c("data.frame", "records", "text")) {

    return.as <- match.arg(return.as)
    if (return.as == "text") {
        return(gene_obj)
    }
    ## Get the records, then call jsonlite:::simplify to convert to a
    ## data.frame
    
    else if (return.as == "data.frame") {
      gene_obj <- .return.as(gene_obj, "records")
      outdf <-jsonlite:::simplify(gene_obj)
      # This expands out any inner columns that may themselves be data frames.
      outdf <- .unnest.df(outdf)
      return(data.frame(DataFrame(outdf, check.names=FALSE), check.names=FALSE))
    }
    ## "DataFrame" is the same as "data.frame" only we convert to the
    ## IRanges class DataFrame and we convert "list" columns to "List"
    ## stances.
#     else if (return.as == "DataFrame") {
#         outdf <- .return.as(gene_obj, "data.frame")
#         return(.df2DF(outdf))
#     }
    else {return(fromJSON(gene_obj, simplifyDataFrame=FALSE))}
}

setGeneric(".request.get", signature=c("mygene"),
            function(mygene, path, params=list()) standardGeneric(".request.get"))

setMethod(".request.get", c(mygene="MyGene"),
            function(mygene, path, params=list()){

    url <- paste(mygene@base.url, path, sep="")
    headers<-c('User-Agent' = sprintf('R-httr_mygene.R/httr.%s', version))
    if (exists('params')){
        if (mygene@debug){
            res <- GET(url, query=params, verbose())
        }
        else{
            res <- GET(url, query=params, config=add_headers(headers))
            }
        }
    if (res$status_code != 200)
        stop("Request returned unexpected status code:\n",
             paste(capture.output(print(res)), collapse="\n"))
    httr::content(res, "text")

})

setGeneric(".request.post", signature=c("mygene"),
            function(mygene, path, params=list()) standardGeneric(".request.post"))

setMethod(".request.post", c(mygene="MyGene"),
            function(mygene, path, params=list()) {

    url <- paste(mygene@base.url, path, sep="")
    headers<-c('Content-Type'= 'application/x-www-form-urlencoded',
            'User-Agent' = sprintf('R-httr_mygene.R/httr.%s', version))
    if (exists('params')){
        if (mygene@debug){
            res <- POST(url, body=params, config=list(add_headers(headers)), verbose())
        }
        else{
            res <- POST(url, body=params, config=list(add_headers(headers)))
            }
        }
    if (res$status_code != 200)
        stop("Request returned unexpected status code:\n",
             paste(capture.output(print(res)), collapse="\n"))
    httr::content(res, "text")
})


.repeated.query <- function(mygene, path, vecparams, params=list(), return.as) {

    verbose <- mygene@verbose
    vecparams.split <- .transpose.nested.list(lapply(vecparams, .splitBySize, maxsize=mygene@step))
    if (length(vecparams.split) <= 1){
        verbose <- FALSE
    }
    vecparams.splitcollapse <- lapply(vecparams.split, lapply, .collapse)
    n <- length(vecparams.splitcollapse)
    reslist <- character(n)
    i <- 1
    repeat {
        if (verbose) {
          message("Querying chunk ", i)
        }
        params.i <- c(params, vecparams.splitcollapse[[i]])
        reslist[[i]] <- .request.post(mygene=mygene, path, params=params.i)
        ## This avoids an extra sleep after the last fragment
        if (i == n){
            break()
        }
        Sys.sleep(mygene@delay)
        i <- i+1
    }
    # This gets the text that would have been returned if we could submit all genes in a single query.
    restext <- .json.batch.collapse(reslist)
    return(restext)
}

setMethod("metadata", c(x="MyGene"), function(x, ...) {
    .return.as(.request.get(x, "/metadata"), "records")
})

setGeneric("getGene", signature=c("mygene"),
            function(geneid, fields = c("symbol","name","taxid","entrezgene"),
            ..., return.as=c("records", "text"), mygene) standardGeneric("getGene"))

setMethod("getGene", c(mygene="MyGene"),
            function(geneid, fields = c("symbol","name","taxid","entrezgene"),
            ..., return.as=c("records", "text"), mygene) {

    return.as <- match.arg(return.as)
    params <- list(...)
    params$fields <- .collapse(fields)
    res <- .request.get(mygene, paste("/gene/", geneid, sep=""), params)
    .return.as(res, return.as=return.as)
})

## If nothing is passed for the mygene argument, just construct a
## default MyGene object and use it.
setMethod("getGene", c(mygene="missing"),
            function(geneid, fields = c("symbol","name","taxid","entrezgene"),
            ..., return.as=c("records", "text"), mygene) {

    mygene <- MyGene()
    callGeneric(geneid, fields, ..., return.as=return.as, mygene=mygene)
})

setGeneric("getGenes", signature=c("mygene"),
            function(geneids, fields = c("symbol","name","taxid","entrezgene"),
            ..., return.as = c("data.frame", "records", "text"), mygene) standardGeneric("getGenes"))

setMethod("getGenes", c(mygene="MyGene"),
            function(geneids, fields = c("symbol","name","taxid","entrezgene"),
            ..., return.as = c("data.frame", "records", "text"), mygene) {

    return.as <- match.arg(return.as)
    if (exists('fields')) {
        params <- list(...)
        params[['fields']] <- .collapse(fields)
        params <- lapply(params, function(x) {str(x);.collapse(x)})
    }
    vecparams <- list(ids=.uncollapse(geneids))
    res <- .repeated.query(mygene, '/gene/', vecparams=vecparams, params=params)
    .return.as(res, return.as=return.as)

})

setMethod("getGenes", c(mygene="missing"),
            function(geneids, fields = c("symbol","name","taxid","entrezgene"),
            ..., return.as = c("data.frame", "records", "text"), mygene) {

    mygene <- MyGene()
    callGeneric(geneids, fields, ..., return.as=return.as, mygene=mygene)
})

setGeneric("query", signature=c("mygene"),
            function(q, ..., return.as=c("data.frame", "records", "text"), mygene) standardGeneric("query"))

setMethod("query", c(mygene="MyGene"),
            function(q, ..., return.as=c("data.frame", "records", "text"), mygene) {

    return.as <- match.arg(return.as)
    params <- list(...)
    params[['q']] <- q
    res <- .request.get(mygene, paste("/query/", sep=""), params)
    if (return.as == "data.frame"){
        return(fromJSON(res))
    }
    else if (return.as == "text"){
        return(.return.as(res, "text"))
    }
    else if (return.as == "records"){
        return(.return.as(res, "records"))
    }

})

setMethod("query", c(mygene="missing"),
            function(q, ..., return.as=c("data.frame", "records", "text"), mygene) {

    mygene <- MyGene()
    callGeneric(q, ..., return.as=return.as, mygene=mygene)
})

setGeneric("queryMany", signature=c("mygene"),
            function(qterms, scopes=NULL, ..., return.as=c("data.frame", 
            "records", "text"), mygene) standardGeneric("queryMany"))

setMethod("queryMany", c(mygene="MyGene"),
            function(qterms, scopes=NULL, ..., return.as=c("data.frame", 
            "records", "text"), mygene){

    return.as <- match.arg(return.as)
    params <- list(...)
    vecparams<-list(q=.uncollapse(qterms))
    if (exists('scopes')){
        params<-lapply(params, .collapse)
        params[['scopes']] <- .collapse(scopes)
         returnall <- .pop(params,'returnall', FALSE)
         params['returnall'] <-NULL
        verbose <- mygene@verbose

        li_query <-c()
        li_missing<-c()
        out <- .repeated.query(mygene, '/query/', vecparams=vecparams, params=params)
        out.li <- .return.as(out, "records")

        for (hits in out.li){
            if (is.null(hits$notfound)){
                li_query<-c(li_query, hits[['query']])
            }
            else if (hits$notfound) {
                li_missing<-c(li_missing, hits[['query']])
            }
        }
        #check duplication hits
        li_cnt<-as.list(table(li_query))
        li_dup<-li_cnt[li_cnt > 1]
        
        if (verbose){
            cat("Finished\n")
            if (length('li_dup')>0){
                sprintf('%f input query terms found dup hits:   %s', length(li_dup), li_dup)
            }
            if (length('li_missing')>0){
                sprintf('%f input query terms found dup hits:   %s', length(li_missing), li_missing)
                }
            }
        out <- .return.as(out, return.as=return.as)
        if (returnall){
            return(list('out'= out, 'dup'=li_dup, 'missing'=li_missing))
        }
        else {
            if (verbose & ((length(li_dup)>=1) | (length(li_missing)>=1))){
                cat('Pass returnall=TRUE to return lists of duplicate or missing query terms.\n')
                }
            return(out)    
        }
    }
})

setMethod("queryMany", c(mygene="missing"),
            function(qterms, scopes=NULL, ...,
            return.as=c("data.frame", "records", "text"), mygene){

    mygene<-MyGene()
    callGeneric(qterms, scopes, ..., return.as=return.as, mygene=mygene)
})

## TODO: Delete below here after migrating everything to S4 ####################################################################

# mygene<-setRefClass("mygene",  fields=c('Url', 'delay', 'step', 'params', 'list','li', 'value', 'geneid', 'geneids', 'fields',
#         'query_fn', 'query_li', 'qr', 'q', 'qterms') , methods =list(

#     initialize=function(Url='http://mygene.info/v2'){
#         .self$Url<<-Url
#         .self$delay<<-1
#         .self$step<<-1000},

#     .pop=function(list, item, default_value=NULL){
#         if (is.null(list[[item]])){
#             return(default_value)}
#         else{
#             value <<- list[[item]]
#             return(value)}},

#     .get=function(Url, params=list()){
#         debug <- .pop(params,'debug', FALSE)
#         params['debug']<<-NULL
#         return_raw <- .pop(params,'return_raw', FALSE)
#         params['return_raw']<<-NULL

#         headers<-c('User-Agent' = sprintf('R-httr_mygene.R/httr.%s', version))
#         if (exists('params')){
#             if (debug){
#                 return(content(GET(Url, query=params, verbose())))}
#             else{
#                 res <- GET(Url, query=params, config=add_headers(headers))
#                 if (res$status_code ==200){
#                     if (return_raw){
#                       return(content(res, 'text'))}
#                     else {
#                       return(fromJSON(content(res, 'text')))}}
#                 else{print(res)}}}},

#     .post=function(Url, params=list()){

#         debug <- .pop(params,'debug', FALSE)
#         params['debug']<<-NULL
#         return_raw <- .pop(params,'return_raw', FALSE)
#         params['return_raw']<<-NULL

#         headers<-c('Content-Type'= 'application/x-www-form-urlencoded',
#             'User-Agent' = sprintf('R-httr_mygene.R/httr.%s', version))
#         if (exists('params')){
#             if (debug){
#                 return(POST(Url, body=params, config=list(add_headers(headers)), verbose()))}
#             else{
#                 res <- POST(Url, body=params, config=list(add_headers(headers)))
#                 if (res$status_code == 200){
#                     if (return_raw){
#                         return(content(res, 'text'))}
#                     else {
#                       out<-content(res, 'text')
#                         return(fromJSON(out, simplifyDataFrame=F))}}
#                 else{print(res)}}}},

#     .format_list=function(list) {
#         out<-paste(list, sep = "", collapse = ",")
#         return(out)},

#     .repeated_query=function(query_fn, query_li, params, verbose=TRUE){
#         result<-list()
#         ql<-length(query_li)
#         if (ql <= .self$step){
#             #No need to do series of batch queries, turn off verbose output
#             verbose = FALSE}

#         for (i in seq(1, ql, by = .self$step)){
#             is_last_loop <- (i+.self$step) >= ql
#             if (verbose){
#                 sprintf('querying %f-%f...', (i+1), min((i+.self$step), ql))}
#             query_result<-query_fn((query_li[i: min(ql, (i+.self$step))]), params)
#             qr<<-append(result, query_result)}
#             if (!is_last_loop & .self$delay){
#                 sys.sleep(.self$delay)}
#             if (verbose){
#                 cat('done.\n')}
#             return(qr)
#             },

#     metadata=function(...){
#         params<<-list(...)
#         .url<-paste(.self$Url, '/metadata', sep = "")
#         return(.self$.get(.url, params))},

#     getgene=function(geneid, fields = 'symbol,name,taxid,entrezgene', ...){
#         # '''Return the gene object for the given geneid.
#         #        This is a wrapper for GET query of "/gene/<geneid>" service.
#         #          @param geneid: entrez/ensembl gene id
#         #          @param fields: fields to return, a list of a comma-sep string
#         #                         if fields=="all", all available fields are returned.
#         #          @param species: optionally, you can pass comma-separated species names
#         #                           or taxonomy ids

#         #        Ref: http://mygene.info/doc/annotation_service.html
#         #     '''
#         if (exists('fields')){
#             params <<- list(...)
#             params[['fields']]<<-fields
#             params<<- lapply(params, function(x) {str(x);.self$.format_list(x)})
#         }
#         .url<-paste(.self$Url, '/gene/', geneid, sep = "")
#         return(.self$.get(.url, params))},

#     .transpose.nested.list=function(li) {
#       ## Assumes that inner names of each element are the same
#       inner.names <- names(li[[1]])
#       setNames(lapply(inner.names, function(i) List(lapply(li, `[[`, i))), inner.names)},

#     .getgenes_inner=function(geneids, params){
#         params[['ids']]<<-.self$.format_list(geneids)
#         .url<-paste(.self$Url, '/gene/', sep = "")
#         return(.self$.post(.url, params))},


#     getgenes=function(geneids, fields = 'symbol,name,taxid,entrezgene', ...){
#         # Return the list of gene object for the given list of geneids.
#      #           This is a wrapper for POST query of "/gene" service.
#      #             @param geneids: a list or comm-sep entrez/ensembl gene ids
#      #             @param fields: fields to return, a list of a comma-sep string
#      #                            if fields=="all", all available fields are returned.
#      #             @param species: optionally, you can pass comma-separated species names
#      #                              or taxonomy ids
#      #             @param filter: alias for fields

#      #          Ref: http://mygene.info/doc/annotation_service.html

#         if (exists('fields')){
#             params <<- list(...)
#             params[['fields']]<<-fields
#             params <<- lapply(params, function(x) {str(x);.self$.format_list(x)})
#         }
#         verbose <- .pop(params,'verbose', TRUE)
#         params['verbose']<<-NULL
#         out<-.self$.repeated_query(.self$.getgenes_inner, geneids, params)
#         return(out)},

#     query=function(q, ...){
#         #         Return  the query result.
#         #         This is a wrapper for GET query of "/query?q=<query>" service.
#         #         @param fields: fields to return, a list of a comma-sep string
#         #                         if fields=="all", all available fields are returned.
#         #         @param species: optionally, you can pass comma-separated species names
#         #                           or taxonomy ids. Default: human,mouse,rat.
#         #         @param size:   the maximum number of results to return (with a cap
#         #                           of 1000 at the moment). Default: 10.
#         #         @param skip:    the number of results to skip. Default: 0.
#         #         @param sort:    Prefix with "-" for descending order, otherwise in ascending order.
#         #                         Default: sort by matching scores in decending order.
#         #         @param entrezonly:  if True, return only matching entrez genes, otherwise, including matching
#         #                              Ensemble-only genes (those have no matching entrez genes).

#         #         Ref: http://mygene.info/doc/query_service.html

#         params<<-list(...)
#         params[['q']] <<- q
#         .url<-paste(.self$Url, '/query/', sep = "")
#         return(.self$.get(.url, params))},

#     .querymany_inner=function(qterms, params){
#         params[['q']]<<-.self$.format_list(qterms)
#         .url<-paste(.self$Url, '/query/', sep = "")
#         return(.self$.post(.url, params))},

#     querymany=function(qterms, scopes=NULL, ...){
#             #     Return the batch query result.
#             #     This is a wrapper for POST query of "/query" service.

#             #     @param qterms: a list of query terms, or a string of comma-separated query terms.
#             #     @param scopes:  type of types of identifiers, either a list or a comma-separated fields to specify type of
#             #                    input qterms, e.g. "entrezgene", "entrezgene,symbol", ["ensemblgene", "symbol"]
#             #                    refer to "http://mygene.info/doc/query_service.html#available_fields" for full list
#             #                    of fields.
#             #     @param fields: fields to return, a list of a comma-sep string
#             #                     if fields=="all", all available fields are returned.
#             #     @param species: optionally, you can pass comma-separated species names
#             #                       or taxonomy ids. Default: human,mouse,rat.
#             #     @param entrezonly:  if True, return only matching entrez genes, otherwise, including matching
#             #                          Ensemble-only genes (those have no matching entrez genes).

#             #     @param returnall:   if True, return a dict of all related data, including dup. and missing qterms
#             #     @param verbose:     if True (default), print out infomation about dup and missing qterms

#             #     Ref: http://mygene.info/doc/query_service.html

#         params<<-list(...)
#         if (exists('scopes')){
#                 params[['scopes']] <<- .self$.format_list(scopes)
#             if ('scope' %in% params){
#                 #allow scope for back-compatibility
#                 params[['scopes']] <<- .self$.format_list(params[['scope']])}
#             if ('fields' %in% params){
#                 params[['fields']] <<- .self$.format_list(params[['fields']])}
#             if ('species' %in% params){
#                 params[['species']] <<- .self$.format_list(params[['species']])}
#             returnall <- .pop(params,'returnall', FALSE)
#             params['returnall']<<-NULL
#             verbose <- .pop(params,'verbose', TRUE)
#             params['verbose']<<-NULL
#             return.as.list <- .pop(params,'return.as.list', FALSE)
#             params['return.as.list']<<-NULL

#             li_missing <-list()
#             li_query <-list()
#             li_cnt <-list()
#             li_dup <-list()
#             out<-.self$.repeated_query(.self$.querymany_inner, qterms, params, verbose=verbose)

#             df <- DataFrame(.self$.transpose.nested.list(out))

#             for (hits in out){
#               if (is.null(hits$notfound)){
#                 li_query<-append(li_query, hits[['query']])}
#               else if (hits$notfound) {
#                 li_missing<-append(li_missing, hits[['query']])}}
#             if (verbose){
#               cat("Finished\n")}
#             #check dup hits
#             li_cnt<-as.list(table(as.character(li_query)))
#             for (hits in li_cnt){
#               if (li_cnt[[hits]] > 1){
#                 li_dup<-append(li_dup, li_cnt[hits])}}

#             if (verbose){
#               if (exists('li_dup')){
#                 sprintf('%f input query terms found dup hits:   %s', length(li_dup), li_dup)}
#               if (exists('li_missing')){
#                 sprintf('%f input query terms found dup hits:   %s', length(li_missing), li_missing)}}
#             if (returnall){
#               return(list('out'= df, 'dup'=li_dup, 'missing'=li_missing))}
#             if (return.as.list){
#               return(out)
#             }
#             else {
#               if (verbose & (exists('li_dup') | exists('li_missing'))){
#                 cat('Pass returnall=TRUE to return lists of duplicate or missing query terms.\n')
#                 return(df)}}}}
#         ))


# mg.getgene<-function(geneid, fields = 'symbol,name,taxid,entrezgene', ...){
#     mg<-mygene$new()
#     return(mg$getgene(geneid=geneid, fields=fields, ...))
# }
# mg.getgenes<-function(geneids, fields = 'symbol,name,taxid,entrezgene', ...){
#     mg<-mygene$new()
#     return(mg$getgenes(geneids=geneids, fields=fields, ...))
# }
# mg.query<-function(q, ...){
#     mg<-mygene$new()
#     return(mg$query(q=q, ...))
# }
# mg.querymany<-function(qterms, scopes=NULL, ...){
#     mg<-mygene$new()
#     return(mg$querymany(qterms=qterms, scopes=scopes, ...))
# }
# mg.metadata<-function(...){
#     mg<-mygene$new()
#     return(mg$metadata(...))}