fillMatrixList  <- function( L, maxLen )
{
    if (is.null(L))     { return(L) }

    # If input is not a list, return original data
    if (!is.list(L))    { return(L) }

    lens <- c()
    for (i in 1:length(L))
    {    
        if (!is.null(L[[i]]))
        {
            lens <- c(lens, getLengths(L[[i]]))
        }
    }

    if (is.null(lens)) { return(NULL) }


    # Set number of columns.
    if (missing(maxLen)) { maxLen <- max(lens, na.rm=TRUE) }


    # If no filling is necessary, return original data
    if (length(unique(lens)) == 1) { return(L) }

    # Enforce maximum length on each subpath of each pathway
    res <- L
    for (i in 1:length(L))
    {
        submat <- L[[i]]
        if (is.null(submat)) { next() }
        replmat <- matrix(0, nrow=nrow(submat), ncol=maxLen)
        for (j in 1:nrow(submat))
        {   
            replmat[j,] <- c(submat[j,], rep(0, maxLen - length(submat[j,])))
        }
        res[[i]] <- replmat
    }

    return(res)
}

unlistToMatrix  <- function( L, mode='rbind' )
{
    if (is.null(L))     { return(L) }

    # If input is already a matrix, return original data.
    if (class(L) == 'matrix') { return(L) }

    # Find what type of data the list holds
    type <- NULL
    for (i in 1:length(L))
    {
        if(!is.null(L[[i]]))
        {
            if (is.vector(L[[i]])) type <- 'vector'
            if (is.matrix(L[[i]])) type <- 'matrix'
            break()
        }
    }

    if (is.null(type)) { return(NULL) }

    # List of vectors of variable size
    if (type == 'vector')
    {
        lens <- sapply(L, function(x) { length(x) })
        M    <- matrix(0, nrow=length(L), ncol=max(lens))
        for (i in 1:length(L)) 
        {
            M[i,1:lens[i]] <- L[[i]] 
        }
        rownames(M) <- names(L)
    }

    # A list of matrices with at least one common dimension size
    # rbind:same number of columns, cbind:same number of rows
    if (type == 'matrix')
    {
        M      <- do.call(mode, L)
        lnames <- names(L) 

        if (length(lnames) > 0)
        {
            lens   <- sapply(L, function(x) 
                                { 
                                    if (!is.null(x)) nrow(x) else 0 
                                } )
            rnames <- vector(mode='numeric', length=nrow(M))
            ctr   <- 0
            for (i in 1:length(lnames))
            {
                if (lens[i] > 0)
                {
                    idx         <- (ctr+1) : (ctr<-ctr+lens[i])
                    rnames[idx] <- rep(lnames[i], lens[i])               
                }
            }
            rownames(M) <- rnames
        }
    }
    return(M)
}

getLengths      <- function( M )
{
    # 
    # Count length of non zero elements in each row of a matrix.
    # Each row must consist of a prefix with non zero elements, 
    # and a suffix of consecutive zeros denote absence of data.
    # Ideal for matrices with a large number of rows.
    #
    esub <- cbind(M, rep(0, nrow(M)))
    esub <- rbind(esub, c(-1, rep(0, ncol(M))))
    df   <- which(diff(which(c(t(esub)) != 0)) - 1 > 0)
    len  <- c(df[1], diff(df))

    return(len)
}


dir.create.rec <- function(fp) 
{
    if( !file.exists(fp) ) 
    {
        dir.create.rec(dirname(fp))
        if ( gsub('\\.', '', fp) == fp  ) 
        {
            dir.create(fp, recursive=FALSE, showWarnings=TRUE) 
        }
    }
} 