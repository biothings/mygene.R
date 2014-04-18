
library(httr)

version <- '0.3'

MyGene<-setRefClass("MyGene",  fields=c('Url', 'delay', 'step', 'params', 'list', 'value', 'geneid', 'geneids', 'fields',
        'query_fn', 'query_li', 'qr', 'q', 'qterms') , methods =list(

    initialize=function(Url='http://mygene.info/v2'){
        .self$Url<<-Url
        .self$delay<<-1
        .self$step<<-1000},

    .pop=function(list, item, default_value=NULL){
        if (is.null(list[[item]])){
            return(default_value)}
        else{
            value <<- list[[item]]
            return(value)}},

    .get=function(Url, params=list()){
    	debug <- .pop(params,'debug', FALSE)
        params['debug']<<-NULL
        return_raw <- .pop(params,'return_raw', FALSE)
        params['return_raw']<<-NULL

        headers<-c('User-Agent' = sprintf('R-httr_mygene.R/httr.%s', version))
        if (exists('params')){
    		res <- GET(Url, query=params, config=add_headers(headers))}
            if (res$status_code ==200){
                if (debug){
                    return(content(GET(Url, query=params, verbose())))}
                if (return_raw){
                    return(content(res, 'text'))}
                else {return(content(res, 'parsed'))}}
            else{print(res)}},

    .post=function(Url, params=list()) {
        debug <- .pop(params,'debug', FALSE)
        params['debug']<<-NULL
        return_raw <- .pop(params,'return_raw', FALSE)
        params['return_raw']<<-NULL

        headers<-c('Content-Type'= 'application/x-www-form-urlencoded',
            'User-Agent' = sprintf('R-httr_mygene.R/httr.%s', version))

        if (exists('params')){

            res <- POST(Url, body=params, config=list(add_headers(headers)))}
            if (res$status_code == 200){
                if (debug){
                    return(content(POST(Url, body=params, config=list(add_headers(headers)), verbose())))}
                if (return_raw){
                    return(content(res, 'text'))}
                else {return(content(res, 'parsed'))}}
            else{print(res)}},

    .format_list=function(list) {
        out<-paste(list, sep = "", collapse = ",")
        return(out)},

    .repeated_query=function(query_fn, query_li, params, verbose=TRUE){
        result<-list()
        ql<-length(query_li)
        if (ql <= .self$step){
            #No need to do series of batch queries, turn off verbose output
            verbose = FALSE}

        for (i in seq(1, ql, by = .self$step)){
            is_last_loop <- (i+.self$step) >= ql
            if (verbose){
                sprintf('querying %f-%f...', (i+1), min((i+.self$step), ql))}
            query_result<-query_fn((query_li[i: min(ql, (i+.self$step))]), params)
            qr<<-append(result, query_result)}
            if (!is_last_loop & .self$delay){
                sys.sleep(.self$delay)}
            if (verbose){
                print('done.')}
            return(qr)
            },

    metadata=function(){
        .url<-paste(.self$Url, '/metadata', sep = "")
        return(.self$.get(.url))},

    getgene=function(geneid, fields = 'symbol,name,taxid,entrezgene', ...){
        # '''Return the gene object for the given geneid.
        #        This is a wrapper for GET query of "/gene/<geneid>" service.
        #          @param geneid: entrez/ensembl gene id
        #          @param fields: fields to return, a list of a comma-sep string
        #                         if fields=="all", all available fields are returned.
        #          @param species: optionally, you can pass comma-separated species names
        #                           or taxonomy ids

        #        Ref: http://mygene.info/doc/annotation_service.html
        #     '''
    	if (exists('fields')){
    		params <<- list(...)
    		params[['fields']]<<-fields
    		params<<- lapply(params, function(x) {str(x);.self$.format_list(x)})
    	}
        .url<-paste(.self$Url, '/gene/', geneid, sep = "")
        return(.self$.get(.url, params))},

    .getgenes_inner=function(geneids, params){
        params[['ids']]<<-.self$.format_list(geneids)
        .url<-paste(.self$Url, '/gene/', sep = "")
        return(.self$.post(.url, params))},


    getgenes=function(geneids, fields = 'symbol,name,taxid,entrezgene', ...){
    	# Return the list of gene object for the given list of geneids.
     #           This is a wrapper for POST query of "/gene" service.
     #             @param geneids: a list or comm-sep entrez/ensembl gene ids
     #             @param fields: fields to return, a list of a comma-sep string
     #                            if fields=="all", all available fields are returned.
     #             @param species: optionally, you can pass comma-separated species names
     #                              or taxonomy ids
     #             @param filter: alias for fields

     #          Ref: http://mygene.info/doc/annotation_service.html

        if (exists('fields')){
    		params <<- list(...)
    		params[['fields']]<<-fields
    		params <<- lapply(params, function(x) {str(x);.self$.format_list(x)})
    	}
        verbose <- .pop(params,'verbose', TRUE)
        params['verbose']<<-NULL
        out<-.self$.repeated_query(.self$.getgenes_inner, geneids, params)
        return(out)},

    query=function(q, ...){
        #         Return  the query result.
        #         This is a wrapper for GET query of "/query?q=<query>" service.
        #         @param fields: fields to return, a list of a comma-sep string
        #                         if fields=="all", all available fields are returned.
        #         @param species: optionally, you can pass comma-separated species names
        #                           or taxonomy ids. Default: human,mouse,rat.
        #         @param size:   the maximum number of results to return (with a cap
        #                           of 1000 at the moment). Default: 10.
        #         @param skip:    the number of results to skip. Default: 0.
        #         @param sort:    Prefix with "-" for descending order, otherwise in ascending order.
        #                         Default: sort by matching scores in decending order.
        #         @param entrezonly:  if True, return only matching entrez genes, otherwise, including matching
        #                              Ensemble-only genes (those have no matching entrez genes).

        #         Ref: http://mygene.info/doc/query_service.html

    	params<<-list(...)
    	params[['q']] <<- q
    	.url<-paste(.self$Url, '/query/', sep = "")
    	return(.self$.get(.url, params))},

    .querymany_inner=function(qterms, params){
        params[['q']]<<-.self$.format_list(qterms)
        .url<-paste(.self$Url, '/query/', sep = "")
        return(.self$.post(.url, params))},

    querymany=function(qterms, scopes=NULL, ...){
            #     Return the batch query result.
            #     This is a wrapper for POST query of "/query" service.

            #     @param qterms: a list of query terms, or a string of comma-separated query terms.
            #     @param scopes:  type of types of identifiers, either a list or a comma-separated fields to specify type of
            #                    input qterms, e.g. "entrezgene", "entrezgene,symbol", ["ensemblgene", "symbol"]
            #                    refer to "http://mygene.info/doc/query_service.html#available_fields" for full list
            #                    of fields.
            #     @param fields: fields to return, a list of a comma-sep string
            #                     if fields=="all", all available fields are returned.
            #     @param species: optionally, you can pass comma-separated species names
            #                       or taxonomy ids. Default: human,mouse,rat.
            #     @param entrezonly:  if True, return only matching entrez genes, otherwise, including matching
            #                          Ensemble-only genes (those have no matching entrez genes).

            #     @param returnall:   if True, return a dict of all related data, including dup. and missing qterms
            #     @param verbose:     if True (default), print out infomation about dup and missing qterms

            #     Ref: http://mygene.info/doc/query_service.html

    	params<<-list(...)
    	if (exists('scopes')){
                params[['scopes']] <<- .self$.format_list(scopes)
            if ('scope' %in% params){
                #allow scope for back-compatibility
                params[['scopes']] <<- .self$.format_list(params[['scope']])}
            if ('fields' %in% params){
                params[['fields']] <<- .self$.format_list(params[['fields']])}
            if ('species' %in% params){
                params[['species']] <<- .self$.format_list(params[['species']])}
        returnall <- .pop(params,'returnall', FALSE)
        params['returnall']<<-NULL
        verbose <- .pop(params,'verbose', TRUE)
        params['verbose']<<-NULL}

        li_missing <-list()
        li_query <-list()
        li_cnt <-list()
        li_dup <-list()
        out<-.self$.repeated_query(.self$.querymany_inner, qterms, params, verbose=verbose)
        for (hits in out){
            if (is.null(hits$notfound)){
                li_query<-append(li_query, hits[['query']])}
            else if (hits$notfound) {
                li_missing<-append(li_missing, hits[['query']])}}

        if (verbose){
            print("Finished")}
        #check dup hits
        li_cnt<-as.list(table(as.character(li_query)))

        for (hits in li_cnt){
            if (li_cnt[[hits]] > 1){
                li_dup<-append(li_dup, li_cnt[hits])}}

        if (verbose){
            if (exists('li_dup')){
                sprintf('%f input query terms found dup hits:   %s', length(li_dup), li_dup)}
            if (exists('li_missing')){
                sprintf('%f input query terms found dup hits:   %s', length(li_missing), li_missing)}}
        if (returnall){
            return(list('out'= out, 'dup'=li_dup, 'missing'=li_missing))}
        else {
            if (verbose & (exists('li_dup') | exists('li_missing'))){
            print('Pass returnall=TRUE to return lists of duplicate or missing query terms.')
            return(out)}}}
        ))


