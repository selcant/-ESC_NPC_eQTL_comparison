# Just making a stand alone script to run GSVA
# I can't seem to get it to work inside the Rmds
# Keep getting this error:
# Error in .mapGeneSetsToFeatures(gset.idx.list, rownames(expr)) : 
# No identifiers in the gene sets could be matched to the identifiers in the expression data.

# Have no idea what is wrong! It was working fine before 09/21/2022

library(here)
load(here("_data","Data_for_GSVA.RData"))
library(GSVA)

gsva_rna <- gsva(  expr = t(expr.npc_rna_upd),
                   gset.idx.list = genesbygo,
                   method ="gsva",
                   kcdf = "none",
                   min.sz = 5, 
                   max.sz = 1000,
                   mx.diff = TRUE
)

gsva_rna_esc_npc <- gsva(  expr = t(shared.expr.merged),
                           genesbygo,
                           method ="gsva",
                           kcdf = "none",
                           min.sz = 5, 
                           max.sz = 1000,
                           mx.diff = TRUE)

save(gsva_rna, gsva_rna_esc_npc, file = here("_data",paste0("Results_from_GSVA_",Sys.Date(),".RData")))

