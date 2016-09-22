
subpathwayTypes <- function(grouping='all')
{
    stream.choices <- c('fwd',
                        'bwd')
    commun.choices <- c('walktrap', 
                        'edge_betweenness', 
                        'fast_greedy', 
                        'leading_eigen', 
                        'infomap', 
                        'louvain')
    compon.choices <- c('decompose', 
                        'max_cliques',
                        'cliques',
                        'coreness')
    topological.choices <- c(
                        'topological.degree', 
                        'topological.betweenness', 
                        'topological.closeness', 
                        'topological.hub_score', 
                        'topological.eccentricity', 
                        'topological.page_rank',
                        'topological.start_nodes')

    functional.choices <- paste0('functional.', .getFunctionalMeasures())
    source.choices <- c(topological.choices, functional.choices)
    
    # Neighborhood
    cases <- expand.grid( stream.choices, 'neighbourhood', source.choices)
    neighborCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})
    # Substream
    cases <- expand.grid( stream.choices, 'stream', source.choices)
    subStreamCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})
    # Linear
    cases <- expand.grid( stream.choices, 'cascade', source.choices)
    linearCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})
    # Community
    cases <- expand.grid( 'community', commun.choices)
    communityCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})
    # Component
    cases <- expand.grid( 'component', compon.choices[1:2])
    compCases_a <- apply(cases, 1, function(x) {paste0(x, collapse='.')})
    cases <- expand.grid( 'component.', paste0(3:9, '-', compon.choices[3]) )
    compCases_b <- apply(cases, 1, function(x) {paste0(x, collapse='')})
    cases <- expand.grid( 'component.', paste0(3:9, '-', compon.choices[4]) )
    compCases_c <- apply(cases, 1, function(x) {paste0(x, collapse='')})
    componentCases <- c(compCases_a, compCases_b, compCases_c)

    supportedMethods <- c(  subStreamCases, 
                            neighborCases, 
                            linearCases, 
                            communityCases, 
                            componentCases)


    if ( is.null(grouping) || grouping == '' ) { return(NULL) }


    if ( length(grouping) == 1 && grepl('^all', grouping) )
    {
        subType <- supportedMethods
        type <- gsub('all.', '', grouping)
        intMeasures <- c('bwd', 'fwd', 'stream', 'neighbourhood', 
                        'cascade', 'community', 'component', 
                        'topological', 'functional', 'DEG')
        extMeasures <- .getExternalMeasures()

        if ( type %in% c(intMeasures, extMeasures) )
        {
            supportedMethods <- subType[grepl(type, subType)]
        }
    }
    if ( length(grouping) == 1 && !grepl('^all', grouping) )
    {
        supportedMethods <- grouping
    }

    return( supportedMethods )
}


.subpathwayAnalysis <- function( edgeList, method=c(), DEgenes, a=3, b=10, 
                                org, verbose=TRUE )
{
    # a: Number of neighbors
    # b: Number of top-nodes 

    supportedMethods <- subpathwayTypes()

    unsupportedOptions <- method[!method %in% supportedMethods]
    if ( length(unsupportedOptions) > 0 )
    {
        message('Option(s) ', unsupportedOptions, ' not supported.')
        out        <- list(NULL)
        names(out) <- paste0('subAnalysis.', method)
        return( out )
    }

    if (verbose)
    {
        message('Performing subpathway analysis ( ', method , ' )...', 
                appendLF = FALSE) 
    }

    # Create an igraph object from an edgelist
    gi <- graph_from_edgelist(as.matrix(edgeList[, 1:2]))

    # Change direction of substream if necessary
    if ( gsub('bwd.', '', method) != method  )
    {
        edgeList <- get.edgelist(gi)
        gi <- graph_from_edgelist(cbind(edgeList[, 2], edgeList[, 1]))
    }
    
    if ( nrow(edgeList) == 0 )
    {
        if (verbose) { message('done.') }
        return (NULL)
    }

    stream.choices <- c('fwd',
                        'bwd')
    commun.choices <- c('walktrap', 
                        'edge_betweenness', 
                        'fast_greedy', 
                        'leading_eigen', 
                        'infomap', 
                        'louvain'
                        )
    compon.choices <- c('decompose', 
                        'max_cliques',
                        'cliques',
                        'coreness'
                        )
    topological.choices <- c(
                        'topological.degree', 
                        'topological.betweenness', 
                        'topological.closeness', 
                        'topological.hub_score', 
                        'topological.eccentricity', 
                        'topological.page_rank',
                        'topological.start_nodes')

    functional.choices <- paste0('functional.', .getFunctionalMeasures())

    source.choices <- c(topological.choices, functional.choices)


    cases <- expand.grid( stream.choices, 'neighbourhood', source.choices)
    neighborCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( stream.choices, 'stream', source.choices)
    subStreamCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( stream.choices, 'cascade', source.choices)
    linearCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( 'community', commun.choices)
    communityCases <- apply(cases, 1, function(x) {paste0(x, collapse='.')})

    cases <- expand.grid( 'component', compon.choices[1:2])
    compCases_a <- apply(cases, 1, function(x) {paste0(x, collapse='.')})
    cases <- expand.grid( 'component.', paste0(3:9, '-', compon.choices[3]) )
    compCases_b <- apply(cases, 1, function(x) {paste0(x, collapse='')})
    cases <- expand.grid( 'component.', paste0(3:9, '-', compon.choices[4]) )
    compCases_c <- apply(cases, 1, function(x) {paste0(x, collapse='')})
    componentCases <- c(compCases_a, compCases_b, compCases_c)


    if ( method %in% neighborCases )
    {
        sourceMeasure <- gsub('fwd.neighbourhood.', '', method)
        sourceMeasure <- gsub('bwd.neighbourhood.', '', sourceMeasure)

        nodes <- .measureToNodes(graph=gi,
                                measure=sourceMeasure,
                                org=org,
                                DEgenes=DEgenes )
        sourceNodes <- nodes[1 : min(b, length(nodes))]
        
        if (length(nodes) == 0)
        { 
            out        <- list(NULL)
            names(out) <- paste0('subAnalysis.', method)
            if (verbose) { message('done.') }
            return(out)
        }

        R <- vector(mode='list', length=length(sourceNodes))
        for (i in seq_len(length(sourceNodes)) )
        {
            R[[i]] <- names(unlist(ego(graph=gi, order=a, 
                                    nodes=sourceNodes[i])))
        }
        names(R) <- paste0('sub', seq_len(length(R)))

        out        <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% subStreamCases )
    {
        sourceMeasure <- gsub('fwd.stream.', '', method)
        sourceMeasure <- gsub('bwd.stream.', '', sourceMeasure)
        nodes <- .measureToNodes(graph=gi,
                                measure=sourceMeasure,
                                org=org,
                                DEgenes=DEgenes )
        sourceNodes <- nodes[1 : min(b, length(nodes))]
        
        if (length(nodes) == 0)
        { 
            out        <- list(NULL)
            names(out) <- paste0('subAnalysis.', method)
            if (verbose) { message('done.') }
            return(out)
        }

        R <- vector(mode='list', length=length(sourceNodes))
        for (i in seq_len(length(sourceNodes)) )
        {
            res <- names(dfs(graph=gi, root=sourceNodes[i], 
                            unreachable=FALSE)[['order']])
            R[[i]] <- res[which(!is.na(res))]
        }
        names(R) <- sourceNodes

        # Keep subpathways with more than two members
        R <- R[which( lapply(R, function(x) { length(x) }) >= 3 )]
        if ( length(R) > 0 ) { names(R) <- paste0('sub', seq_len(length(R)) ) }

        out        <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% linearCases )
    {
        sourceMeasure <- gsub('fwd.cascade.', '', method)
        sourceMeasure <- gsub('bwd.cascade.', '', sourceMeasure)

        # Find all source nodes
        nodes <- .measureToNodes(graph=gi,
                                measure=sourceMeasure,
                                org=org,
                                DEgenes=DEgenes )
        sourceNodes <- nodes[1 : min(b, length(nodes))]

        if (length(nodes) == 0)
        { 
            out        <- list(NULL)
            names(out) <- paste0('subAnalysis.', method)
            if (verbose) { message('done.') }
            return(out)
        }

        # Simplify graph
        adjmat <- get.adjacency(gi)
        gi <- graph.adjacency(triu(adjmat))

        # Find all destination nodes
        destinNodes <- names(sort(which(degree(gi, mode='out', 
                                loops=FALSE) == 0), decreasing=TRUE))
        sourceNodes <- sourceNodes[!sourceNodes %in% destinNodes]
        destinNodes <- destinNodes[!destinNodes %in% sourceNodes]
                
        # Find all linear subpathways between each pair of start/end nodes
        N <- length(sourceNodes)*length(destinNodes)
        lpaths <- vector(mode='list', length=N)
        lnames <- vector(mode='numeric', length=N)
        cnt <- 1

        lpaths <- NULL
        if ( N > 0 )
        {
            for ( sourceNode in sourceNodes )
            {
                for ( destinNode in destinNodes )
                {
                    lpaths[[cnt]] <- all_simple_paths( gi, 
                                                from = sourceNode, 
                                                to = destinNode,
                                                mode="out")
                    lnames[cnt] <- paste(sourceNode, destinNode, sep='-')
                    cnt <- cnt + 1
                }
            }

            combinations <- expand.grid( sourceNodes, destinNodes )
            names(lpaths) <- lnames
            names(lpaths) <- paste0('sub', seq_len(length(lpaths)) )

            # Remove NULL results and flatten first level of results 
            lpaths <- do.call(c, lpaths)
            # Keep subpathways with more than two members
            lpaths <- lpaths[which(lapply(lpaths, length) > 2)]
            # Keep genes ids
            lpaths <- lapply(lpaths, function(x) { names(x) } )
        }

        out        <- list(lpaths)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% communityCases )
    {
        if ( method == 'community.walktrap' )
        { 
            gr <- cluster_walktrap(gi)
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) { R[[i]] <- (gr[[i]]) }
            }
        }
        if ( method == 'community.edge_betweenness' )
        { 
            gr <- cluster_edge_betweenness(gi)
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) { R[[i]] <- (gr[[i]]) }
            }
        }
        if ( method == 'community.fast_greedy' )
        { 
            gr <- cluster_fast_greedy(as.undirected(gi))
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) { R[[i]] <- (gr[[i]]) }
            }
        }
        if ( method == 'community.leading_eigen' )
        { 
            gr <- cluster_leading_eigen(as.undirected(gi))
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) { R[[i]] <- (gr[[i]]) }
            }
        }
        if ( method == 'community.infomap' )
        { 
            gr <- cluster_infomap(as.undirected(gi))
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) { R[[i]] <- (gr[[i]]) }
            }
        }
        if ( method == 'community.louvain' )
        { 
            gr <- cluster_louvain(as.undirected(gi))
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) { R[[i]] <- (gr[[i]]) }
            }       
        }

        if ( length(R) > 0 ) { names(R) <- paste0('sub', seq_len(length(R)) ) }

        out <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    if ( method %in% componentCases )
    {
        if ( method == 'component.max_cliques' )
        { 
            gr <- max_cliques(as.undirected(gi), a, b)
            R <- vector( mode='list', length=length(gr) )
            if ( length(gr) > 0 )
            {
                for ( i in seq_len(length(gr)) ) 
                                { R[[i]] <- names(unlist(gr[i])) }
            }
        }

        if ( method == 'component.decompose')
        {
            gr  <- decompose(gi)
            R   <- lapply(gr, function(x) { as_data_frame(x, what="edges") } )
        }

        if ( gsub('-cliques', '', method) !=  method )
        {
            p <- gsub('-cliques', '', method)
            k <- as.numeric(gsub('component.', '', p))

            g <- as_graphnel( gi )
            ksubs <- kCliques(ugraph(g))[paste0(k, '-cliques')][[1]]
            ksubs <- lapply(ksubs, function(x) { matrix(x, nrow=1) } )

            # Extract valid pairs
            edgeList <- as_data_frame(gi, what="edges")
            
            # Consider only actual interactions between genes.
            if ( length(ksubs) > 0 )
            {
                ksubs <- .unlistToMatrix(.fillMatrixList(ksubs))
                R <- vector(mode='list', length=nrow(ksubs))
                for (j in seq_len(nrow(ksubs)) )
                {
                    # Edge list indexes
                    r1  <- as.numeric(is.element( edgeList[,1], ksubs[j, ]) )
                    r2  <- as.numeric(is.element( edgeList[,2], ksubs[j, ]) )
                    idx <- which((r1 + r2) == 2)
                    if ( length(idx) > 0 )
                    { 
                        R[[j]] <- unique(as.vector(t(edgeList[idx, ])))
                    }
                }
            }else{ R <- NULL}

        }
        if ( gsub('-coreness', '', method) !=  method )
        {
            p <- gsub('-coreness', '', method)
            p <- as.numeric(gsub('component.', '', p))

            corDist  <- coreness(gi)
            R <- names(corDist)[which(corDist == p)]
        }
        if ( length(R) > 0 ) { names(R) <- paste0('sub', seq_len(length(R)) ) }


        out <- list(R)
        names(out) <- paste0('subAnalysis.', method)
    }
    
    if (verbose) { message('done.') }


    return(out)
}


.getFunctionalNodes <- function( graph, targets, org )
{
    # Extract nodes from graph and keep term related ones
    graphGenes <- names(V(graph))
    # Unique target genes in entrez 
    uGenesFromTargets <- unique(as.vector(targets))
    uGenesFromTargets <- unname(.changeAnnotation(data=uGenesFromTargets, 
                                org='hsa', choice='HGNCtoEntrez'))
    uGenesFromTargets <- uGenesFromTargets[!is.na(uGenesFromTargets)]
    # Keep targets intersecting with graph genes
    targetGenesInGraph <- graphGenes[graphGenes %in% uGenesFromTargets]
    targetGenesInGraph <- unname(.changeAnnotation(data=targetGenesInGraph, 
                                org='hsa', choice='entrezToHGNC'))

    # Filter target genes with graph genes
    idx <- matrix(targets %in% targetGenesInGraph, nrow=nrow(targets))*1
    targets[!idx] <- NA

    nodes <- names(sort(table(targets), decreasing=TRUE))
    # Change to entrez gene annotation
    nodes <- unname(.changeAnnotation(data=nodes, org=org, 
                    choice='HGNCtoEntrez'))

    return (nodes)
}


.measureToNodes <- function ( graph, measure, org, DEgenes=NULL )
{
    if ( measure == 'topological.degree' ) 
    { 
        nodes <- names(sort(degree(graph), decreasing=TRUE))
    }
    if ( measure == 'topological.betweenness' ) 
    { 
        nodes <- names(sort(betweenness(graph), decreasing=TRUE))
    }
    if ( measure == 'topological.closeness' ) 
    {
        nodes <- names(sort(closeness(graph), decreasing=TRUE))
    }
    if ( measure == 'topological.hub_score' ) 
    { 
        nodes <- names(sort(hub_score(graph)[['vector']], decreasing=TRUE))
    }
    if ( measure == 'topological.eccentricity' ) 
    { 
        nodes <- names(sort(eccentricity(graph), decreasing=TRUE))
    }
    if ( measure == 'topological.page_rank' ) 
    { 
        nodes <- names(sort(page_rank(graph)[['vector']], decreasing=TRUE))
    }
    if ( measure == 'topological.start_nodes' ) 
    {
        func <- function(x) 
                { which(degree(x, mode='in', loops=FALSE) == 0) } 
        nodes <- names(sort(func(graph), decreasing=TRUE))
    }
    if ( measure == 'functional.DEG' ) 
    { 
        nodes <- names(sort(DEgenes, decreasing=FALSE))
    }

    if ( measure %in% paste0('functional.', .getExternalMeasures() ) )
    {
        if ( org != 'hsa' )
        { 
            out <- list(NULL)
            names(out) <- paste0('subAnalysis.', measure)
            return( out )
        }

        # Find genes with most occurences in the selected term
        measure <- gsub('functional.', '', measure)
        targets <- .loadTermData( type=measure )
        nodes <- .getFunctionalNodes(graph=graph, targets=targets, org=org)
    }

    return( nodes )
}



#
# Data-related functions
#
.changeAnnotation <- function(data, org, choice)
{
    if ( org == 'hsa' )
    {
        dir <- system.file('extdata//Data', package='DEsubs')

        load(paste(dir, 'libraryEntrezToHGNC.RData', sep='//'), 
            e <- new.env())
        libraryEntrezToHGNC <- e[['libraryEntrezToHGNC']]
        if ( choice == 'entrezToHGNC' )
        {
            data <- libraryEntrezToHGNC[data]    
        }
        if ( choice == 'HGNCtoEntrez' )
        {
            libraryHGNCtoEntrez <- names(libraryEntrezToHGNC)
            names(libraryHGNCtoEntrez) <- libraryEntrezToHGNC
            data <- libraryHGNCtoEntrez[data]    
        }
    }

    return(data)
}

.getExternalMeasures <- function()
{
    defaultReferences <- c( 'KEGG',
                            'GO_bp',
                            'GO_cc',
                            'GO_mf',
                            'Disease_OMIM',
                            'Disease_GAD',
                            'Drug_DrugBank',
                            'miRNA',
                            'TF')

    # Search for new gene sets
    files <- list.files(cache[['datDir']])
    otherFiles <- c('libraryEntrezToHGNC.RData', 
                    'edgeLists.RData',
                    'libraryEntrezToExternalNomenclature.RData')
    files <- files[-which(files %in% otherFiles)]

    # Find if any of the default files have been deleted
    idx <- which(!paste0(defaultReferences, '.RData') %in% files)
    if ( length(idx) > 0 )
    {
        message( 'References ', paste0(files[idx], collapse=', '), 
                    ' missing.')
    }
    supportedReferences <- gsub('.RData', '', files)

    return( supportedReferences )
}

.getFunctionalMeasures <- function()
{

    return( c('DEG', .getExternalMeasures()) )
}


#
# Base data
#

.loadTermData <- function( type )
{
    dir <- cache[['datDir']]
    supportedTerms <- .getExternalMeasures()

    if ( type %in% .getExternalMeasures() )
    {
        iFile <- paste0( dir, '//' ,type, '.RData' )
        load(iFile, e <- new.env())        
        targetsPerClass <- e[['targetsPerClass']]
    }else{
        message('Type ', type, ' not supported.')
    }

    if ( is.null(targetsPerClass) )
        { message('Set ', type, ' does not contain valid targers.') }

    return( targetsPerClass )
} 



