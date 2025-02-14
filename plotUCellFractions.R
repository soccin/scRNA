require(tidyverse)

argv=yaml::read_yaml("pass_02b_PARAMS.yaml")

THETA=0.1
UCELLFILE="uCellScores_lung_Gardner_v1_.csv"

md=read_csv(UCELLFILE) %>%
    select(CellID,SampleID,matches("_UC$")) %>%
    rename_all(~gsub("_UC$","",.))

ct1=md %>% mutate(Type=ifelse(Basal.Cells>THETA,"Basal.Cells",NA))
ct2=md %>%
    select(CellID,3:ncol(.)) %>%
    select(-Basal.Cells) %>%
    gather(Type,UCell,-CellID) %>%
    group_by(CellID) %>%
    slice_max(UCell) %>%
    ungroup %>%
    filter(UCell>THETA)

ct=ct1 %>%
    left_join(ct2,by="CellID") %>%
    mutate(Type=case_when(!is.na(Type.x) ~ Type.x, !is.na(Type.y) ~ Type.y, T ~ "unknown")) %>%
    select(-Type.x,-Type.y)

pg=ct %>%
    ggplot(aes(SampleID,fill=Type)) +
        theme_light(20) +
        geom_bar(position="fill") +
         xlab(NULL) +
         ylab(NULL) +
         scale_y_continuous(labels = scales::percent_format(accuracy = 1))

mTag=basename(UCELLFILE) %>% tools::file_path_sans_ext() %>% gsub("_$","",.) %>% gsub("uCell[^_]*_","",.) %>% gsub("_"," ",.)

p2=pg + scale_fill_brewer(palette="Paired") + labs(title=mTag,subtitle=argv$PROJNAME)

pFile=cc(argv$PROJNAME,"moduleFractionsUCell",mTag,".pdf") %>% gsub(" ","_",.)
pdf(file=pFile,height=8.5,width=11)
print(p2)
dev.off()

pctType=ct %>%
    count(SampleID,Type) %>%
    group_by(SampleID) %>%
    mutate(PCT=n/sum(n)) %>%
    ungroup %>%
    select(-n) %>%
    spread(SampleID,PCT,fill=0)

tFile=cc(argv$PROJNAME,"moduleFractionsUCell",mTag,".xlsx") %>% gsub(" ","_",.)
openxlsx::write.xlsx(pctType,tFile)

