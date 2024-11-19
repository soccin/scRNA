require(tidyverse)
require(readxl)

getCellTypeGeneLists<-function(listOfAllGenes) {

    typeTable=read_xlsx("markerGenesV2.xlsx") %>%
        select(Type,Genes) %>%
        separate_rows(Genes,sep=",") %>%
        mutate(Genes=gsub(" ","",Genes) %>% str_to_title(.)) %>%
        filter(Genes %in% listOfAllGenes) %>%
        group_split(Type)

    types=list()
    for(ii in seq(len(typeTable))) {
        types[typeTable[[ii]][[1]][1]]=typeTable[[ii]][2]
    }

    types

}
