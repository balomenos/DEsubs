cache <- new.env()


.onLoad <- function(libname, pkgname)
{

    #
}


.onAttach <- function(libname, pkgname)
{
    .buildDirectories()
    #
}


.onUnload <- function(libname) 
{

    # 
}


.buildDirectories  <- function()
{
    # Set the default location of the base directory
    path <- switch(.Platform[['OS.type']], unix = path.expand("~"),
                    windows= file.path(gsub("\\\\", "/",
                    Sys.getenv("USERPROFILE")), "AppData"))
    opt <- getOption("DEsubs_CACHE", paste0(path, '/DEsubs'))
    baseDir <- Sys.getenv("DEsubs_CACHE", opt)

    usrDir <- paste0(baseDir, '//User')
    outDir <- paste0(baseDir, '//Output')
    datDir <- paste0(baseDir, '//Data')

    dir.create(baseDir, showWarnings=FALSE, recursive=TRUE)
    dir.create(outDir, showWarnings=FALSE, recursive=TRUE)

    cache[['baseDir']] <- baseDir
    cache[['usrDir']] <- usrDir
    cache[['outDir']] <- outDir
    cache[['datDir']] <- datDir

    # Copy demo files from package directory to user directort
    file.copy(  from=system.file('extdata//Data', package='DEsubs'), 
                to=file.path(baseDir), 
                overwrite = FALSE, recursive = TRUE, copy.mode = TRUE)

    file.copy(  from=system.file('extdata//User', package='DEsubs'), 
                to=file.path(baseDir), 
                overwrite = FALSE, recursive = TRUE, copy.mode = TRUE)

    return(invisible())
}

