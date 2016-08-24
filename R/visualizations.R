
# Gene level visualizations
matrixVisualization <- function( data, type, title='', colors, width, height, 
                                outfile, export )
{
    # Specify base dir for default storage
    dir <- cache$outDir 
    # Create user specified directory structure 
    if ( !missing(outfile) && (outfile != '') && ('pdf' %in% export) ) 
    { 
        dir.create.rec(outfile) 
    }

    if ( type == 'heatmap')
    {
        if (missing(colors))
            { colors <- colorRampPalette(c("white", "red"))(100) }

        if ( missing(outfile) ) 
        { 
            outfile <- paste0(dir, '//heatmap_', title, '.pdf') 
        }

        if ( 'pdf' %in% tolower(export) && (!'plot' %in% tolower(export) ) )
        { 
            pdf(outfile, width=width, height=height, onefile=FALSE)
        }

        hmap <- pheatmap(data, 
                cluster_rows=FALSE, 
                cluster_cols=FALSE,
                cellwidth=40,
                cellheight=15,
                legend=TRUE,
                color=colors)
        if ( 'plot' %in% tolower(export) && 'pdf' %in% tolower(export) )
        {         
            dev.copy2pdf(file=outfile, width=width, height=height)
        }
    }
    if ( type == 'barplot')
    {
        if (missing(colors))
            { colors <- c('#C7EDFCFF')}

        if ( 'plot' %in% tolower(export) )
        { 
            par(mar=c(5, 1, 0, 1))
        }
        
        if ( missing(outfile) ) 
            { outfile <- paste0(dir, '//barplot', title, '.pdf') }
        
        if ( 'pdf' %in% tolower(export) && (!'plot' %in% tolower(export) ) )
        { 
            pdf(outfile, width=width, height=height, onefile=FALSE) 
        }
        midpoints <- barplot(data, 
                            col = colors[1], 
                            horiz=TRUE, 
                            xlab='-log10(Q.value)',
                            ylab='Genes',
                            axes = TRUE, axisnames = FALSE, 
                            border=NA)
        text(x=0, y=midpoints-0.1, pos=4, cex=0.75, labels=names(data))

        if ( 'plot' %in% tolower(export) && 'pdf' %in% tolower(export) )
        {         
            dev.copy2pdf(file=outfile, width=width, height=height)
        }
    }
    if ( type == 'dotplot')
    {
        if ( missing(outfile) || outfile == '' ) 
            { outfile <- paste0(dir, '//dotplot_', title[2], '.pdf') }
        
        if ( 'pdf' %in% tolower(export) && (!'plot' %in% tolower(export) ) )
        { 
            pdf(outfile, width=width, height=height, onefile=FALSE) 
        }
        colors <- colorRampPalette(c("white", "red"))(100)
        print(qplot(data[, 1], data[, 2], data = data, colour = data[, 3], 
                size=I(3))
                + theme_classic() + geom_point()
                + scale_colour_gradientn(colours=colors)
                + labs(x=title[1], y=title[2], colour=title[3])
                + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, 
                                                    hjust=1))
                )
        if ( 'plot' %in% tolower(export) && 'pdf' %in% tolower(export) )
        {         
            dev.copy2pdf(file=outfile, width=width, height=height)
        }
    }

    if ( 'pdf' %in% tolower(export) && (!'plot' %in% tolower(export) ) )
    { 
        dev.off() 
    }
}


# Subpathway to graph visualization
subpathwayToGraph <- function(data, submethod, subname, colors, size, 
                                export, width, height, outfile, verbose)
{

    if ( missing(export) ) { export <- 'pdf' }
    if ( missing(verbose) ) { verbose <- TRUE }
    if ( missing(size) ) { size <- c(4,4) }

    if ( verbose )
        { message('Visualizing the subpathway as a graph...', 
            appendLF = FALSE) }

    doPlot <- ( length(export) == 1 ) && ( 'plot' %in% tolower(export) )
    if ( !missing(outfile) && !doPlot ) 
    {
        if ( gsub('\\*', '', outfile) != outfile ) 
            { outfile <- gsub('\\*', 'pdf', outfile) }
        dir.create.rec(outfile)
    }

    # Assign node colors according to pValue of node
    if (missing(colors))
    {
        # Assign node colors according to pValue of node
        colors <- c( '#FFFFFFFF', '#FFF0F0FF', '#FF0000FF', '#FF0000FF', 
                '#FF5C5CFF', '#FF5C5CFF', '#FFC9C9FF', '#FFC9C9FF',
                '#FF8585FF', '#FF8585FF')
        # colors <- colorRampPalette(c("white", "red"))(100) 
    }

    exportTypes <- c('pdf', 'plot',  'edgelist', 'json', 'gml', 
                        'ncol', 'lgl', 'graphml','dot')
    idx <- which(!export %in% exportTypes)
    if ( length(idx) > 0)
    {
        message('Types ', paste(export[idx], collapse=', '), ' not supported.')
        message('Availiable types are ', 
                                    paste(exportTypes, collapse=', '), '.')
    }

    edgeList <- data$edgeList
    supportedMethods <- subpathwayTypes()

    unsupportedOptions <- submethod[!submethod %in% supportedMethods]
    if ( length(unsupportedOptions) > 0 )
    {
        message('Option(s) ', unsupportedOptions, ' not supported.')
    }
    if ( length(submethod) > 1 )
    {
        message('Please select a single subpathway type for visualization.')
    }

    submethodName <- paste0('subAnalysis.', submethod)
    subpathway.nodelist <- data[[submethodName]][[subname]]

    if ( is.null(subpathway.nodelist) )
    {
        message('Subpathway ', subname, ' does not exist.')
        return(NULL)
    }

    # Subpathway (nodelist) to edgelist
    idx1 <- which(edgeList[, 1] %in% subpathway.nodelist)
    idx2 <- which(edgeList[, 2] %in% subpathway.nodelist)
    idx <- intersect(idx1, idx2)
    subpathway.edgelist <- edgeList[idx, ]

    DEgenes<- data$DEgenes
    subpathway.nodelist <- DEgenes[names(DEgenes) %in% subpathway.nodelist]
    subpathway.nodelist <- data.frame(  'genes'=names(subpathway.nodelist),
                                        'degree'=unname(subpathway.nodelist),
                                        stringsAsFactors=FALSE )
    idx <- order(as.numeric(subpathway.nodelist[, 1]))
    subpathway.nodelist <- subpathway.nodelist[idx, ]
    org <- data$org

    if (org == 'hsa')
    {
        # Change to gene names
        subpathway.edgelist[, 1] <- changeAnnotation(   
                                    data=subpathway.edgelist[, 1], 
                                    org=org, choice='entrezToHGNC' )
        subpathway.edgelist[, 2] <- changeAnnotation(   
                                    data=subpathway.edgelist[, 2], 
                                    org=org, choice='entrezToHGNC' )
    }

    # Create the network
    net <- graph.data.frame(subpathway.edgelist, directed=TRUE)
    net <- simplify(net, remove.multiple = FALSE, remove.loops = TRUE)  
    scores <- order(subpathway.nodelist$'degree')
    scores <- floor(scores/max(scores)*10)/10
    subpathway.nodelist.norm <- cbind(subpathway.nodelist, 'score'=scores)
    idx <- order(as.vector(subpathway.nodelist.norm[, 3]))
    subpathway.nodelist.norm <- subpathway.nodelist.norm[idx, ]
    idx <- as.numeric( cut(subpathway.nodelist.norm[, 3], 6))

    finalColors <- colors[idx]
    V(net)$color <- finalColors

    # Assign weight according to correlation
    edgeWeights <- vector(mode='numeric', length=nrow(subpathway.edgelist))
    edgeWeights <- rep(1, length(edgeWeights))
    shuffleIdx <- order(as.vector(abs(subpathway.edgelist[, 5])))
    subpathway.edgelist <- subpathway.edgelist[shuffleIdx, ]
    reverseMapper <- shuffleIdx
    names(reverseMapper) <- 1:length(shuffleIdx)
    reverseMapper <- as.numeric(names(sort(reverseMapper)))
    idx <- as.numeric( cut(abs(subpathway.edgelist[, 5]), 6))
    edgeWeights <- edgeWeights + idx*1.5
    subpathway.edgelist.scored <- cbind(subpathway.edgelist, 
                                'weight'=rep(1,nrow(subpathway.edgelist)))


    idx_pos <- which(subpathway.edgelist[, 5] >= 0.7)
    idx_neg <- which(subpathway.edgelist[, 5] <= -0.7)
    idx <- unique(c(idx_neg, idx_pos))
    subpathway.edgelist.scored[idx, 5] <- edgeWeights[idx]


    # Assign edge color according to correlation type (positive or negative)
    edgeColors <- vector(mode='numeric', length=nrow(subpathway.edgelist))
    edgeColors <- rep('black', length(edgeColors))
    edgeColors[idx_pos] <- 'darkgreen'
    edgeColors[idx_neg] <- 'red'

    # Restore to the order corresponding to the graph
    subpathway.edgelist.scored <- subpathway.edgelist.scored[reverseMapper, ]
    edgeColors <- edgeColors[reverseMapper]

    E(net)$arrow.size <- .5
    E(net)$color <- edgeColors
    E(net)$width <- subpathway.edgelist.scored[, 6]

    if ( 'plot' %in% tolower(export) )
    {
        set.seed(2345)
        par(mar=c(0, 0, 0, 1))
        plot(net, 
            layout=layout.davidson.harel, 
            vertex.label.dist=1,
            vertex.label.degree=pi/2)
    }

    if ( 'pdf' %in% tolower(export) )
    {
        dir <- cache$outDir 

        if ( missing(outfile) )
            { out <- paste0(dir, '//subplot_', subname, '.pdf')  } 
        if ( !missing(outfile) )
            { out <- outfile } 

        pdf(out, width=width, height=height, onefile=TRUE)
        set.seed(2345)

        plot(net, 
            layout=layout.davidson.harel, 
            vertex.label.dist=1,
            vertex.label.degree=pi/2) 

        dev.off()
    }

    rownames(subpathway.edgelist) <- 1:nrow(subpathway.edgelist)


    #
    # Other exports
    #
    dir <- cache$outDir
    subpathway.edgelist[, 5] <- floor((subpathway.edgelist[, 5]*100))/100
    subpathway.edgelist <- subpathway.edgelist[, c('E1', 'E2', 'corr')]
    g <- graph.data.frame(subpathway.edgelist, directed=TRUE)

    if ( 'json' %in% tolower(export) )
    {
        if ( missing(outfile) )
            { out <- paste0(dir, '//subpath_', subname, '.json') }
        if ( !missing(outfile) )
            { out <- outfile }

        if ( file_ext(out) != 'json' )
        {
            out <- gsub(file_ext(out), 'json', out)
        }
        subjson <- toJSON( subpathway.edgelist)
        write(file=out, subjson)
    }
    if ( 'edgelist' %in% tolower(export) )
    {
        if ( missing(outfile) )
            { out <- paste0(dir, '//subpath_', subname, '.tsv') }
        if ( !missing(outfile) )
            { out <- outfile }

        if ( file_ext(out) != 'tsv' )
        {
            out <- gsub(file_ext(out), 'tsv', out)

        }
        write.table(file=out, subpathway.edgelist, quote=FALSE, 
                row.names=FALSE, col.names=TRUE)
    }

    otherExportTypes <- c( 'ncol', 'lgl', 'graphml', 'gml', 'dot')
    fileTypes <- c('.ncol', '.lgl', '.xml', '.gml', '.dot')
    names(fileTypes) <- otherExportTypes    

    for (type in export)
    {
        if ( !tolower(type) %in% otherExportTypes ) { next() }
        if ( missing(outfile) )
            { out <- paste0(dir, '//subpath_', subname, fileTypes[type]) }
        if ( !missing(outfile) )
            { out <- outfile }

        if ( file_ext(out) != fileTypes[type] )
        { out <- gsub(file_ext(out), gsub('\\.', '', fileTypes[type]) , out) }

        subgml <- write.graph( g, format=type, file=out)
    }

    if ( verbose )
        { message('done', appendLF = TRUE) }

    return ( invisible() )
}

doCirclize <- function(mat, colors, outfile, agap, bgap, cgap, degree, 
                        a, b, export)
{

    # doCirclize

    if ( 'plot' %in% export )
    { 
        par(mar=c(0, 0, 0, 0))
    }
    
    if ( 'pdf' %in% export && (!'plot' %in% export ) )
    { 
        pdf(outfile)
    }


    order <- c(rownames(mat), rev(colnames(mat)))

    if (missing (agap))   { agap <- 20 }    
    if (missing (bgap))   { bgap <- 20 }    
    if (missing (cgap))   { cgap <- 50 }
    if (missing (degree)) { degree <- 230 }

    if ( agap == 0 )
    {
        gd <- c(rep(2, nrow(mat)-1), cgap,
                rep(2, ncol(mat)-1), cgap)
    }
    if ( agap > 0 && bgap == 0 )
    {
        k1 <- floor((nrow(mat)-1)/2) 
        k2 <- ceiling((nrow(mat)-1)/2)

        if (k1 == 0 || k2 == 0)
        {
            gd <- c(rep(2, nrow(mat)-1), cgap,
                    rep(2, ncol(mat)-1), cgap)
        }
        if (k1 > 0 && k2 > 0)
        {
            gd <- c(rep(2, k1), agap, 
                    rep(2, k2-1), cgap, 
                    rep(2, ncol(mat)-1), cgap)
        }
    }
    if ( agap > 0 && bgap > 0 )
    {
        k1 <- floor((nrow(mat)-1)/3) 
        k2 <- ceiling((nrow(mat)-1)/3)
        k3 <- floor((nrow(mat)-1)/3) 
        gd <- c(rep(2, k1-1), agap, 
                rep(2, k2), bgap, 
                rep(2, k3-1), cgap, 
                rep(2, ncol(mat)-1), cgap)
    }

    if ( missing(colors) )
    {
        colors <- rainbow(nrow(mat))
    }

    gcol <- NULL
    gcol[rownames(mat)] <- colors
    gcol[colnames(mat)] <- 'lightgray'
    par(mfrow=c(1,1), mar=c(0, 0, 0, 0))
    
    tryCatch(
        circos.par(start.degree = degree, gap.degree = gd,
                canvas.xlim=c(-a, b), canvas.ylim=c(-a, b))
    ,error=
        circos.par(start.degree = degree, 
            gap.degree = c(rep(2, nrow(mat)-1), cgap, 
                            rep(2, ncol(mat)-1), cgap),
            canvas.xlim=c(-a, b), canvas.ylim=c(-a, b))
    )
    
    chordDiagram(mat,order=order, grid.col = gcol, transparency = 0.5, 
    annotationTrack="grid", annotationTrackHeight=0.01, preAllocateTracks=1)

    # Since default text facing in `chordDiagram` is fixed, 
    # we need to manually add text in track 1
    for(si in get.all.sector.index()) 
    {
        xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 1)
        ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 1)

        circos.text(
            mean(xlim), ylim[1], si, 
            facing = "clockwise", 
            adj = c(0, 0),
            niceFacing = TRUE, 
            cex = 0.7, 
            col = "black", 
            sector.index = si, 
            track.index = 1
        )
    }

    circos.clear()  

    if ( 'plot' %in% tolower(export) && 'pdf' %in% tolower(export) )
    {         
        dev.copy2pdf(file=outfile)
    }

    if ( 'pdf' %in% tolower(export) && (!'plot' %in% tolower(export) ) )
    { 
        dev.off() 
    }
}


