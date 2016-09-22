.fillMatrixList  <- function( data, maxLen )
{
    # Zero pads a list of matrices (data) each having different number of 
    # columns. If the maximum number of columns is availiable, providing 
    # is as an argument speeds up computation time (maxlen).
    if (is.null(data))     { return(data) }

    # If input is not a list, return original data
    if (!is.list(data))    { return(data) }

    lens <- c()
    for (i in seq_len(length(data)) )
    {    
        if (!is.null(data[[i]]))
        {
            lens <- c(lens, .getLengths(data[[i]]))
        }
    }

    if (is.null(lens)) { return(NULL) }


    # Set number of columns.
    if (missing(maxLen)) { maxLen <- max(lens, na.rm=TRUE) }


    # If no filling is necessary, return original data
    if (length(unique(lens)) == 1) { return(data) }

    # Enforce maximum length on each subpath of each pathway
    res <- data
    for (i in seq_len(length(data)) )
    {
        submat <- data[[i]]
        if (is.null(submat)) { next() }
        replmat <- matrix(0, nrow=nrow(submat), ncol=maxLen)
        for (j in seq_len(nrow(submat)) )
        {   
            replmat[j,] <- c(submat[j,], rep(0, maxLen - length(submat[j,])))
        }
        res[[i]] <- replmat
    }

    return(res)
}

.unlistToMatrix  <- function( data, mode='rbind' )
{
    # Reshapes a list of matrices to one matrix in a full vectorized fashion
    if (is.null(data))     { return(data) }

    # If input is already a matrix, return original data.
    if (class(data) == 'matrix') { return(data) }

    # Find what type of data the list holds
    type <- NULL
    for (i in seq_len(length(data)) )
    {
        if(!is.null(data[[i]]))
        {
            if (is.vector(data[[i]])) type <- 'vector'
            if (is.matrix(data[[i]])) type <- 'matrix'
            break()
        }
    }

    if (is.null(type)) { return(NULL) }

    # List of vectors of variable size
    if (type == 'vector')
    {
        lens <- sapply(data, function(x) { length(x) })
        M    <- matrix(0, nrow=length(data), ncol=max(lens))
        for (i in seq_len(length(data)) ) 
        {
            M[i,1:lens[i]] <- data[[i]] 
        }
        rownames(M) <- names(data)
    }

    # A list of matrices with at least one common dimension size
    # rbind:same number of columns, cbind:same number of rows
    if (type == 'matrix')
    {
        M      <- do.call(mode, data)
        lnames <- names(data) 

        if (length(lnames) > 0)
        {
            lens   <- sapply(data, function(x) 
                                { 
                                    if (!is.null(x)) nrow(x) else 0 
                                } )
            rnames <- vector(mode='numeric', length=nrow(M))
            ctr   <- 0
            for (i in seq_len(length(lnames)) )
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

.getLengths      <- function( data )
{
    # 
    # Count length of non zero elements in each row of a matrix.
    # Each row must consist of a prefix with non zero elements, 
    # and a suffix of consecutive zeros denote absence of data.
    # Ideal for matrices with a large number of rows.
    #
    esub <- cbind(data, rep(0, nrow(data)))
    esub <- rbind(esub, c(-1, rep(0, ncol(data))))
    df   <- which(diff(which(c(t(esub)) != 0)) - 1 > 0)
    len  <- c(df[1], diff(df))

    return(len)
}


.dir.create.rec <- function(file) 
{
    if( !file.exists(file) ) 
    {
        .dir.create.rec(dirname(file))
        if ( gsub('\\.', '', file) == file  ) 
        {
            dir.create(file, recursive=FALSE, showWarnings=TRUE) 
        }
    }
} 