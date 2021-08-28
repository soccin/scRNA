require(biomaRt)
require(Seurat)
require(tidyverse)
require(yaml)

cc.genes=map(cc.genes.updated.2019,as_tibble) %>% bind_rows(.id="Group") %>% rename(Human=value)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

genesV2 = getLDS(
    attributes = c("hgnc_symbol"),
    filters = "hgnc_symbol",
    values = cc.genes$Human,
    mart = human,
    attributesL = c("mgi_symbol"),
    martL = mouse,
    uniqueRows=T)

cc.genes=left_join(cc.genes,genesV2,by=c(Human="HGNC.symbol")) %>% rename(Mouse.biomaRt=MGI.symbol)

human.gene=read_csv("hoy/Jax_Ensembl_Biomart_mouse_to_human_gene_20190524.csv")

db.order <- c("Jax/Ensembl/Biomart",
                "Jax/Biomart", "Jax/Ensembl", "Ensembl/Biomart",
                "Jax", "Enembl", "Biomart")

gene.list=cc.genes
gene.list <- gene.list %>% mutate(Gene = Human)

res=select(gene.list, Gene=Human) %>%
    left_join(human.gene, by=c("Gene"="human")) %>%
    mutate(DB = factor(DB, levels=db.order)) %>%
    filter(!is.na(mouse)) %>%
    group_by(Gene) %>%
    arrange(DB) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(gene.list, ., by="Gene")

cc.genes.mouse.v2=split(res,res$Group) %>% map(.,"mouse")

write_yaml(cc.genes.mouse.v2,"cc.gene.mouse.v2.yaml")

