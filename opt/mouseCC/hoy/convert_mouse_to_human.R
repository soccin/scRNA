convert_human_to_mouse <- function(gene.list, colname="Gene"){
  # add this `Gene` colname for easier left_join
  gene.list <- gene.list %>% mutate(Gene := !!sym(colname))
  
  human.gene <- read_csv("Jax_Ensembl_Biomart_mouse_to_human_gene_20190524.csv")
  db.order <- c("Jax/Ensembl/Biomart", 
                "Jax/Biomart", "Jax/Ensembl", "Ensembl/Biomart", 
                "Jax", "Enembl", "Biomart")
  
  res <- dplyr::select(gene.list, Gene = !!sym(colname)) %>% 
    left_join(., human.gene, by=c("Gene"="human")) %>% 
    mutate(DB = factor(DB, levels=db.order)) %>%
    filter(!is.na(mouse)) %>%
    group_by(Gene) %>% 
    arrange(DB) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    left_join(gene.list, ., by="Gene") %>%
    mutate(mouse = ifelse(is.na(mouse), first_UP(Gene), mouse)) %>% 
    dplyr::select(-Gene, -DB) %>%  
    dplyr::select(Gene = mouse, everything()) 
  return(res)
}

convert_mouse_to_human <- function(gene.list, colname="Gene"){
  # add this `Gene` colname for easier left_join
  gene.list <- gene.list %>% mutate(Gene := !!sym(colname))
  
  human.gene <- read_csv("Jax_Ensembl_Biomart_mouse_to_human_gene_20190524.csv")
  db.order <- c("Jax/Ensembl/Biomart", 
                "Jax/Biomart", "Jax/Ensembl", "Ensembl/Biomart", 
                "Jax", "Enembl", "Biomart")
  
  res <- dplyr::select(gene.list, Gene = !!sym(colname)) %>% 
    left_join(., human.gene, by=c("Gene"="mouse")) %>% 
    mutate(DB = factor(DB, levels=db.order)) %>%
    filter(!is.na(human)) %>%
    group_by(Gene) %>% 
    arrange(DB) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    left_join(gene.list, ., by="Gene") %>%
    mutate(human = ifelse(is.na(human), toupper(Gene), human)) %>% 
    dplyr::select(-Gene, -DB) %>%  
    dplyr::select(Gene = human, everything()) 
  return(res)
}


