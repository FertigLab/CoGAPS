read.data <- function(file, format=c("jags","bugs"))
{
    format <- match.arg(format)
    switch(format, "jags"=read.jagsdata(file), "bugs"=read.bugsdata(file))
}

read.jagsdata <- function(file)
{
  e <- new.env()
  eval(parse(file), e)
  return(as.list(e))
}

read.bugsdata <- function(file)
{
    bugs.dat <- dget(file)
    for (n in names(bugs.dat)) {
        if (!is.null(dim(bugs.dat[[n]]))) {
            dim(bugs.dat[[n]]) <- rev(dim(bugs.dat[[n]]))
            bugs.dat[[n]] <- aperm(bugs.dat[[n]])
        }
    }
    return(bugs.dat)
}
    

