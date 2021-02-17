
git.describe<-function(){

    suppressPackageStartupMessages(require(git2r))

    # system("git describe --always --dirty")

    mostRecentTagIdx=rev(order(unlist(lapply(git2r::tags(),function(x){x$tagger$when$time}))))[1]
    mostRecentTag=git2r::tags()[[mostRecentTagIdx]]
    commits=git2r::commits()

    paste(
        mostRecentTag$name,
        which(lapply(commits,function(x){x$sha})==mostRecentTag$target),
        paste0("g",substr(commits[[1]]$sha,1,7)),
        sep="-"
        )

}