#!/usr/bin/env Rscript
# Running on Sumner
# singularity exec docker://rocker/tidyverse:4.1.0 R
# Rscript to get LOD scores at peaks from scan files 
# Then merge those values across 2 datasets: ESC RNA, NPC RNA

library(magrittr)
library(tidyverse) # load tiyverse packages
library(assertthat) # for assert_that
library(qtl2)
library(GenomicRanges)
select<-dplyr::select
rename<-dplyr::rename
source("/projects/munger-lab/projects/DO_mNPC/rna_seq/scripts/qtl2geno_splitMap.R")

## ESC RNA data
load("/projects/munger-lab/projects/DO_mESC/rna_seq/qtl_mapping/total_gene_expression/eqtl_grid69k_pe/DO_mESC_paired_eQTL_peaks.RData")
peaks.esc_rna <- peaks %>% 
  select( -ci_lo, -ci_hi)
rm(peaks)

## NPC RNA data
load("/projects/munger-lab/projects/DO_mNPC/rna_seq/qtl_mapping/total_gene_expression/eqtl_grid69k_pe/DO_mNPC_paired_eQTL_peaks.RData")
peaks.npc_rna <- peaks %>% 
  select( -ci_lo, -ci_hi)
rm(peaks)

# load gene annotations
load("/projects/munger-lab/projects/DO_mNPC/rna_seq/ENSMUSGid_to_symbol_v91.RData") # gene.info
all.genes <- gene.info 

### Annotating peaks
# Let's add physical pos of peaks
# Again, I modified the code chunk below from Dan
chroms <- c(as.character(1:19), "X")
interp_bp <- function(df) {
  df <- arrange(df, peak_chr, peak_cM)
  peak_gpos <- select(df, peak_chr, peak_cM)
  chr <- peak_gpos$peak_chr
  f <- factor(chr, chroms)
  peak_gcoord_list <- split(peak_gpos$peak_cM, f)
  peak_pcoord_list <- qtl2::interp_map(peak_gcoord_list, gmap, pmap)
  df$interp_bp_peak <- unsplit(peak_pcoord_list, f)
  df
}
# Making a map with basepair positions -- the pos column is in cM not in bp!!
# pmap - has bp attached to markers instead of cM
map_dat2 <- map_dat %>%
  separate(marker, into=c('chrom', 'pos_bp'), convert=T, remove=F) %>%
  mutate(n=1:n()) %>% as_tibble()
pmap <- split_map(dplyr::select(map_dat2, marker,
                                chr, pos_bp) %>% as.data.frame() %>%
                    tibble::remove_rownames() %>%
                    tibble::column_to_rownames('marker'))

peaks.esc_rna.wbp <- interp_bp(peaks.esc_rna)  %>%
  mutate( interp_bp_peak = ifelse( interp_bp_peak ==3e+06, (3e+06)+1, interp_bp_peak))
peaks.npc_rna.wbp <- interp_bp(peaks.npc_rna)  %>%
  mutate( interp_bp_peak = ifelse( interp_bp_peak ==3e+06, (3e+06)+1, interp_bp_peak))


# I am matching by chromosome + phenotype (gene name used in qtl mapping)
# That let's me find the overlap easier
# Then, I identify peaks that are further apart than 5Mb + I will seperate the NPC peaks and add as rows to be included in the matching later
# esc_rna + npc_rna
merged.peaks.esc.npc <- full_join(peaks.esc_rna.wbp, peaks.npc_rna.wbp, 
                                  by=c("phenotype","peak_chr"), 
                                  suffix=c(".esc_rna",".npc_rna")) %>%
  mutate(bp_diff = abs(interp_bp_peak.esc_rna-interp_bp_peak.npc_rna)) %>% # there are QTL for the same gene on the same chr but further >5Mb, which should be added as a separate peak/row 
  mutate(match=ifelse(!is.na(bp_diff), TRUE,NA )) %>% 
  mutate(match=ifelse( match ==TRUE & bp_diff >10e06, FALSE,match)) %>%
  mutate(match=ifelse(is.na(match)&!is.na(lod.esc_rna),"esc_rna",match )) %>%
  mutate(match=ifelse(is.na(match)&!is.na(lod.npc_rna),"npc_rna",match ))
# Fix the match == FALSE ones! 
# add the npc/esc columns as new rows seperately + leave the other as NA + add the lod scores from scans. 
merged.peaks.esc.npc.fixed <- filter(merged.peaks.esc.npc, match !=FALSE) %>% 
  rbind(., (filter(merged.peaks.esc.npc, match==FALSE) %>%
              mutate_at( vars(ends_with(".esc_rna")),function(x){ x=as.double(NA)} )
  ))%>%
  rbind(., (filter(merged.peaks.esc.npc, match==FALSE) %>%
              mutate_at( vars(ends_with(".npc_rna")),function(x){ x=as.double(NA)} )
  ))

# check dims to make sure!
assertthat::assert_that(dim(filter(merged.peaks.esc.npc, match !=FALSE))[1]+
                          dim(filter(merged.peaks.esc.npc.fixed, match ==FALSE & is.na(lod.esc_rna)))[1]+
                          dim(filter(merged.peaks.esc.npc.fixed, match ==FALSE & is.na(lod.npc_rna)))[1]
                        == dim(merged.peaks.esc.npc.fixed)[1])


# Function for adding missing lod scores
add_lods <- function(peaks, scans1, scans2, r=2.5e06, map=map_dat2){
  #peaks : merged peaks data frame, FALSE matches added as new rows. 
  #scans1: ESC_rna 
  #scans2: NPC_rna 
  #r: defines the region to look around the peak position for overlap
  
  for(i in 1:dim(peaks)[1]){
    phenotype <- peaks$phenotype[i] #this is the gene name - will be column name in scans
    
    if( peaks$match[i]==TRUE ){ # if both lods are there move on to the next gene
      #print( paste(i,"common"))
      next
    }
    if( peaks$match[i]==FALSE ){
      if(!is.na(peaks$lod.esc_rna[i]) & is.na(peaks$lod.npc_rna[i]) & phenotype %in% colnames(scans2) ){ # if ESC score is there but NPC is missing fill in the NPC
        
        peak_chr <- peaks$peak_chr[i]
        peak_bp  <- peaks$interp_bp_peak.esc_rna[i]
        region   <- c(peak_bp-r, peak_bp+r)
        region.map <- map %>% filter(., chr==peak_chr & pos_bp > region[1] & pos_bp < region[2])
        region.lod <- scans2[region.map$marker , phenotype]
        max.lod <- max(region.lod)
        max.marker <- names(region.lod)[which(region.lod==max.lod)]
        if(length(max.marker >1)){max.marker<-max.marker[1]}
        max.map <- filter(map, chr==peak_chr & marker == max.marker) # get the matchin pos_bp from the map
        
        # Now, let's fill in the gaps in the peaks data frame
        # change values in peaks
        peaks$lod.npc_rna[i] <- max.lod
        peaks$peak_cM.npc_rna[i] <- max.map$pos
        peaks$interp_bp_peak.npc_rna[i] <- max.map$pos_bp
        
      }
      if(is.na(peaks$lod.esc_rna[i]) & !is.na(peaks$lod.npc_rna[i]) & phenotype %in% colnames(scans1)){ # if NPC score is there but ESC is missing fill in the ESC
        
        peak_chr <- peaks$peak_chr[i]
        peak_bp  <- peaks$interp_bp_peak.npc_rna[i]
        region   <- c(peak_bp-r, peak_bp+r)
        region.map <- map %>% filter(., chr==peak_chr & pos_bp > region[1] & pos_bp < region[2])
        region.lod <- scans1[region.map$marker , phenotype]
        max.lod <- max(region.lod,na.rm=T)
        max.marker <- names(region.lod)[which(region.lod==max.lod)]
        if(length(max.marker >1)){max.marker<-max.marker[1]}
        max.map <- filter(map, chr==peak_chr & marker == max.marker) # get the matchin pos_bp from the map
        
        # Now, let's fill in the gaps in the peaks data frame
        # change values in peaks
        peaks$lod.esc_rna[i] <- max.lod
        peaks$peak_cM.esc_rna[i] <- max.map$pos
        peaks$interp_bp_peak.esc_rna[i] <- max.map$pos_bp
      }
    }
    
  }
  return(peaks)
}
#################################################################################################################################
#################################################################################################################################
############################ Filling in missing lod scores ######################################################################
# load scans.
load("/projects/munger-lab/projects/DO_mNPC/rna_seq/qtl_mapping/total_gene_expression/eqtl_grid69k_pe/DO186_mNPC_paired_eQTL_scans.RData")
load("/projects/munger-lab/projects/DO_mESC/rna_seq/qtl_mapping/total_gene_expression/eqtl_grid69k_pe/DO185_mESC_paired_eQTL_scans.RData")


message("filling in the LOD scores")

merged.peaks.esc.npc.filled <- add_lods(merged.peaks.esc.npc.fixed,
                                        scans1=esc.rna.scans,
                                        scans2= npc.scans,
                                        r=5e06)

# check that the fill is working!
assertthat::assert_that( dim(filter(merged.peaks.esc.npc.fixed, is.na(lod.esc_rna)))[1] >
                           dim(filter(merged.peaks.esc.npc.filled, is.na(lod.esc_rna)))[1]
)
assertthat::assert_that( dim(filter(merged.peaks.esc.npc.fixed, is.na(lod.npc_rna)))[1] >
                           dim(filter(merged.peaks.esc.npc.filled, is.na(lod.npc_rna)))[1]
)

#########################################################
# annotating peaks
# let's change the match == false to the correct match
peaks.esc.npc.rna.wgeneinfo <- merged.peaks.esc.npc.filled %>%  
  # update using case_when
  mutate(
    match = case_when(
      match == "esc_rna" ~ "esc_rna",
      match == "npc_rna" ~ "npc_rna",
      match == TRUE ~ "shared",
      match == FALSE & (lod.esc_rna >  lod.npc_rna | is.na(lod.npc_rna))  ~ "esc_rna",
      match == FALSE & (lod.esc_rna <  lod.npc_rna | is.na(lod.esc_rna))   ~ "npc_rna"
    )
  ) %>% 
  dplyr::rename(ensembl_gene_id=phenotype) %>%
  left_join(.,all.genes) %>%
  dplyr::mutate(gene_chrom=chrom, gene_start=start, gene_end=end, gene_strand=strand)

# Let's add local - distant
# Look at local and distant peaks
# esc_rna + npc_rna
peaks.esc.npc.rna <- peaks.esc.npc.rna.wgeneinfo %>% 
  mutate(same_chrom=peak_chr == gene_chrom) %>%
  mutate(midpoint=(gene_start + gene_end)/2) %>%
  mutate(pos_within_10Mb.esc_rna=abs(interp_bp_peak.esc_rna - midpoint) < 10e6,
         pos_within_10Mb.npc_rna=abs(interp_bp_peak.npc_rna - midpoint) < 10e6)
peaks.esc.npc.rna$local.esc_rna <- with(peaks.esc.npc.rna, same_chrom & pos_within_10Mb.esc_rna)
peaks.esc.npc.rna$local.npc_rna <- with(peaks.esc.npc.rna, same_chrom & pos_within_10Mb.npc_rna)

# prep some stuff for plotting:
uchr <- c(as.character(1:19), "X")
cl <- dplyr::select(map_dat2, chr, pos_bp) %>% group_by(chr) %>%
  summarize(len=max(pos_bp))
clp <- with(cl, setNames(len, chr))
chrom_lens <- setNames(as.numeric(clp[uchr]), uchr)
chrom_lens_offset <- cumsum(chrom_lens) - chrom_lens
chrom_lens_midpt <- chrom_lens_offset + chrom_lens/2

# esc_rna + npc_rna
peaks.esc.npc.rna  <- peaks.esc.npc.rna  %>% 
  mutate(cumsum_bp.esc_rna = interp_bp_peak.esc_rna + chrom_lens[peak_chr],
         cumsum_bp.npc_rna = interp_bp_peak.npc_rna + chrom_lens[peak_chr])
peaks.esc.npc.rna$cumsum_bp_peak.esc_rna <- peaks.esc.npc.rna$interp_bp_peak.esc_rna+chrom_lens_offset[peaks.esc.npc.rna$peak_chr]
peaks.esc.npc.rna$cumsum_bp_peak.npc_rna <- peaks.esc.npc.rna$interp_bp_peak.npc_rna+ chrom_lens_offset[peaks.esc.npc.rna$peak_chr]
peaks.esc.npc.rna$cumsum_bp_gene <- peaks.esc.npc.rna$midpoint + chrom_lens_offset[peaks.esc.npc.rna$gene_chrom]
peaks.esc.npc.rna  <- peaks.esc.npc.rna  %>% arrange(peak_chr)


# let's save the comparison tables. I can use them locally to summarize results in an Rnotebook. 
table.note <-paste0("This data has been prepared by SA using 01_Compare_DO_Peaks_onlyNPC.R using the interactive que on the cluster on ", Sys.Date())
#save(table.note, peaks.esc.npc.rna, file=paste0("/projects/munger-lab/projects/DO_eQTL_comparison/data/", "peaks_comparison_10Mb_ESC_NPC.RData"))
save(table.note, peaks.esc.npc.rna, file=paste0("/projects/munger-lab/projects/DO_mNPC/ESC_NPC_eQTL_comparison/_data/", "peaks_comparison_10Mb_ESC_NPC.RData"))
