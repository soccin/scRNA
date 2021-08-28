
git.describe<-function(REPO="."){

    system2(
        "git",
        c(
            paste0("--git-dir=",file.path(REPO,".git")),
            paste0("--work-tree=",REPO),
            "describe","--tags", "--always", "--long", "--dirty='-UNCOMMITED'"),
        stdout=T
        )

    # suppressPackageStartupMessages(require(git2r))

    # system("git describe --always --dirty")

    # mostRecentTag=git2r::tags()[[1]]
    # mostRecentTagName=names(git2r::tags())[1]
    # commits=git2r::commits()

    # paste(
    #     mostRecentTagName,
    #     which(lapply(commits,function(x){x$sha})==mostRecentTag$sha),
    #     paste0("g",substr(commits[[1]]$sha,1,7)),
    #     sep="-"
    #     )

}