
params=yaml::read_yaml("pass_02b_PARAMS.yaml")

THETA=0.1
argv=commandArgs(trailing=T)
if(len(argv)!=1) {
    cat("\n\n   plotUCellFractions uCellMetaDataFile.csv\n\n")
    quit()
}
UCELLFILE=argv[1]

require(tidyverse)

md0=read_csv(UCELLFILE)
md=md0 %>%
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

p2=pg + scale_fill_brewer(palette="Paired") + labs(title=mTag,subtitle=params$PROJNAME)

pFile=cc(params$PROJNAME,"moduleFractionsUCell",mTag,".pdf") %>% gsub(" ","_",.)
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

tFile=cc(params$PROJNAME,"moduleFractionsUCell",mTag,".xlsx") %>% gsub(" ","_",.)
openxlsx::write.xlsx(pctType,tFile)


sc=md0 %>%
    select(CellID,SampleID,SCT_snn_res.0.5,matches("_UC$")) %>%
    gather(Module,Score,matches("_UC")) %>%
    mutate(Module=gsub("_UC$","",Module))

scoreBySample=sc %>%
    group_by(SampleID,Module) %>%
    summarize(Score=mean(Score)) %>%
    spread(Module,Score)

scoreByCluster=sc %>%
    group_by(SCT_snn_res.0.5,Module) %>%
    summarize(Score=mean(Score)) %>%
    spread(Module,Score)

tFile=cc(params$PROJNAME,"moduleMeansUCell",mTag,".xlsx") %>% gsub(" ","_",.)
openxlsx::write.xlsx(list(ScoreByCluster=scoreByCluster,ScoreBySample=scoreBySample),tFile)

