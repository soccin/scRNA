require(tidyverse)
argv=commandArgs(trailing=T)
oo=readRDS(argv[1])
oo@meta.data %>%
    rownames_to_column("CellID") %>%
    tibble %>%
    write_csv("metadata.csv")
