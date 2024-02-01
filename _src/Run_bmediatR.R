# running bmediatR on distant eQTL

# load libraries
# for data handling + plotting
library(tidyverse)
library(ggpubr)
# for qtl mapping+ annotation + mediation
library(qtl2)
library(GenomicRanges)
library(intermediate)
library(bmediatR)
# re-assign some functions explicitly
select <- dplyr::select 
rename <- dplyr::rename 
summarise <- dplyr::summarise

# load data
# mapping objects
load("/projects/munger-lab/projects/DO_mNPC/ESC_NPC_eQTL_comparison/_data/DO_mNPC_data_for_bmediatr_2023-03-22.RData")

# get distant eQTL
npc_eqtl_dist <- peaks.npc_rna.wEffs %>% 
  filter( lod.npc_rna > 7.5, local.npc_rna ==F)
# get distant eQTL with a good mediator ESC or NPC
npc_eqtl_dist_top_meds <- eqtl_npc_rna_meds %>% 
  mutate( type = "npc") %>% 
  rbind( eqtl_esc_rna_meds %>% mutate( type ="esc")) %>% 
  inner_join(
    npc_eqtl_dist %>% select( target.id = ensembl_gene_id, qtl.chr = peak_chr, target.lod = lod.npc_rna)
  ) %>% 
  filter( mediation.lod < 0.5*(target.lod),
          !str_detect(target.symbol, "-ps"),
          mediator.symbol != "Cwc22") %>% # get high lod drop examples and filter pseudogenes
  group_by( target.id, target.lod, qtl.chr) %>% 
  slice_min( mediation.lod, n = 2) # getting top 2 best mediators to test with bmediatR
  

# run bmediatr
# create an empty list to store results
results <- c()

#create a for loop to run mediation with `bmediatR`
run_bmediatr<- function(expr, qtl, med_expr, med_annot){
  gene_expr
  
  
}
  gene_expression <- exprZ.npc_rna[,data$target.id[i]]
  marker <- map_dat2 %>%
    mutate(diff = abs(pos_bp - data$interp_bp_peak.npc_rna[i])) %>%
    slice_min(diff)
  mediator <- exprZ.esc_rna[,colnames(exprZ.esc_rna) !=data$target.id[i]]
  genotype <- pull_genoprobpos(
    genoprobs = probs.npc_rna,
    marker = marker$marker
  )
  bmediatr_scan <- bmediatR(
    y = gene_expression,
    M = mediator,
    X = genotype, 
    Z = covar.npc_rna
  )
  results[[i]] <- bmediatr_scan


save(results, file = here("data","bmediatr_results_mNPC_eQTL_wESC.RData"))


# save

