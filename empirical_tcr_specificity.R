# ==== set up ====
packages <- c("Seurat", "magrittr", "dplyr", "reshape2", "pheatmap",
              "ggplot2", "ggrepel", "ggsci", "ggpubr",
              "RColorBrewer", "wesanderson", "viridis",
              "gridExtra", "grid", "cowplot")

invisible(lapply(packages, library, character.only = TRUE))

source("~/oak/scripts/sc_utils.R")

# ==== load data ====
cfi <- readRDS(file="r_objects/cfi.rds")
tissue_tet <- readRDS(file="tissue_tetramer/r_objects/tissue_tet.rds")
load(file="r_objects/all_color_pals.RData")

ps_pal <- c(brewer.pal(n = 3, "Reds")[c(3,2)], "#9fcb92", "grey60")
names(ps_pal) <- levels(cfi$ps_category)

# ==== map tissue multimer+ dataset onto total CD3+ dataset ====
DefaultAssay(cfi) <- "integrated" # ref
DefaultAssay(tissue_tet) <- "RNA" # query

anchors <- FindTransferAnchors(reference=cfi, query=tissue_tet, 
                                dims=1:30, reference.reduction="pca")

tissue_tet <- MapQuery(anchorset = anchors, reference = cfi, query = tissue_tet,
                       refdata = list(pheno_cluster = "pheno_cluster"), 
                       reference.reduction = "pca", reduction.model = "umap")

tissue_tet$predicted.pheno_cluster <- factor(tissue_tet$predicted.pheno_cluster,
                                            levels = levels(cfi$pheno_cluster))

DimPlot(tissue_tet, reduction = "ref.umap", group.by = "predicted.pheno_cluster", 
        label = TRUE, cols = cluster_pal,
        label.size = 3, repel = TRUE, raster=T) + 
    NoLegend() + ggtitle("UMAP projection")

# cluster annotation
tt_cluster_annotation <- c(`0`="TRM", `1`="GZMK_high", `2`="Memory", 
                           `3`="Ex", `4`="HSP_high", `5`="Effector", 
                           `6`="MT_high", `7`="HSP_high2", `8`="Prolif_Ex",
                           `9`="CD4_contam", `10`="CC_high", `11`="Ex_ISG", 
                           `12`="MAIT", `13`="Treg_contam",
                           `14`="KIR", `15`="Effector_ZEB2"
                          )

tt_cluster_levels <- c("Effector", "Effector_ZEB2", "Memory", "TRM", 
                       "GZMK_high", "HSP_high", "HSP_high2", "CC_high",
                       "Ex", "Ex_ISG", "Prolif_Ex", "MT_high", "KIR", "MAIT", 
                       "CD4_contam", "Treg_contam"
                      )

tissue_tet$pheno_cluster <- factor(tt_cluster_annotation[tissue_tet$seurat_clusters],
                                   levels = tt_cluster_levels)

tt_cluster_pal <- colorRampPalette(brewer.pal(9, "Set1"))(16)[c(4,5,3,11,10,2,6,13,8,9,1,12,7,16,15,14)]
names(tt_cluster_pal) <- tt_cluster_levels

DimPlot(tissue_tet, group.by="pheno_cluster", reduction = "umap", 
        label=T, cols = tt_cluster_pal, raster = T)

# filter out CD4s (small amount of contamination from sort)
Idents(tissue_tet) <- "pheno_cluster"
tissue_tet <- subset(tissue_tet, idents = c("Treg_contam", "CD4_contam"), invert=T)

DefaultAssay(tissue_tet) <- "integrated"

tissue_tet <- ScaleData(tissue_tet, verbose = FALSE) %>% 
    RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)

DimPlot(tissue_tet, group.by="pheno_cluster", reduction = "umap", 
        label=T, cols = tt_cluster_pal, raster = T)

ggplot(subset(tissue_tet[[]], !pheno_cluster %in% c("CD4_contam", "Treg_contam")), 
              aes(x=patient)) +
    geom_bar(aes(fill=predicted.pheno_cluster), position="fill") +
    scale_fill_manual(values=cluster_pal, name="Projected cluster", limits=force) +
    coord_cartesian(ylim=c(0,1), expand=F) +
    labs(x="", y="Proportion of cells") +
    theme_bw()

# MSK1263 only
tissue_tet_1263 <- subset(tissue_tet, patient=="MSK1263")

DefaultAssay(tissue_tet_1263) <- "integrated"

tissue_tet_1263 <- ScaleData(tissue_tet_1263, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>% 
    RunUMAP(dims = 1:30)

# Fig S10C: tissue multimer cluster UMAP
DimPlot(tissue_tet_1263, group.by="pheno_cluster", reduction = "umap", label=T, 
        cols = tt_cluster_pal, raster=T)
ggsave("plots/tissue_tet_cluster_umap_msk1263.pdf", width=8, height=6)

# Fig S10D: projected vs. multimer+ cluster annotations
confusion <- table(tissue_tet_1263[[c("predicted.pheno_cluster", "pheno_cluster")]])
confusion <- confusion[rowSums(confusion)!=0,colSums(confusion)!=0]

pheatmap::pheatmap(confusion[c("MAIT", "CD8-EFF", "CD8-PROLIF-EXH", "CD8-EXH", "CD8-GZMK", "CD8-TRM"),], 
                   main = "Cluster projection", angle=45,
                   cluster_rows = F,
                   scale = "column",
                   file = "plots/tissue_tet_cfi_mapquery_confusion_horizontal_msk1363.pdf", 
                   height=3.2, width=6,
                   treeheight_row = 15, treeheight_col = 10)

# Fig S10E: projected cluster bar plot
ggplot(subset(tissue_tet_1263[[]], !pheno_cluster %in% c("CD4_contam", "Treg_contam")), 
              aes(x=patient)) +
    geom_bar(aes(fill=predicted.pheno_cluster), position="fill") +
    scale_fill_manual(values=cluster_pal, name="Projected cluster", limits=force) +
    coord_cartesian(ylim=c(0,1), expand=F) +
    labs(x="", y="Proportion of cells") +
    theme_bw()

ggsave("plots/tissue_tet_projected_phenotypes_barplot_msk1263.pdf", width=3, height=4)

# ==== read empirical neopeptide-specific TCRs ====
# PB traditional MANAFEST bulk TCR
pb_t_manafest_clones <- readRDS(file="r_objects/pb_t_manafest_clones.rds")
pb_t_manafest_ps_trbs <- pb_t_manafest_clones[pb_t_manafest_clones$specificity == "NeoAg", "pat_TRB_aa"]
pb_t_manafest_vs_trbs <- pb_t_manafest_clones[pb_t_manafest_clones$specificity == "ViralAg", "pat_TRB_aa"]

# PB modified MANAFEST scTCR
pb_m_manafest_clones_df <- readRDS(file="blood_modified_manafest/r_objects/pb_m_manafest_tetpos_clones_df.rds")
pb_m_manafest_ps_trbs <- unique(pb_m_manafest_clones_df$pat_TRB_aa)

# tissue dextramer bulk TCR
bulk_tissue_dex_ps_trbs <- readRDS(file = "r_objects/na_clones_bulk_tissue_dex.rds")

# tissue tetramer scTCR
sc_tissue_tet_df <- readRDS(file="tissue_tetramer/r_objects/tissue_tet_ps_clones.rds")
sc_tissue_tet_ps_trbs <- unique(sc_tissue_tet_df$pat_TRB_aa)

# remove ViralAg specific clones from traditional MANAFEST
bulk_tissue_dex_ps_trbs_novs <- setdiff(bulk_tissue_dex_ps_trbs, pb_t_manafest_vs_trbs)
sc_tissue_tet_ps_trbs_novs <- setdiff(sc_tissue_tet_ps_trbs, pb_t_manafest_vs_trbs)

# ==== define empirical TCR specificity of clones =====
# all peptide-specific specific clones
# remove ViralAg specific clones from traditional MANAFEST
all_ps_trbs <- setdiff(unique(c(sc_tissue_tet_ps_trbs, 
                                bulk_tissue_dex_ps_trbs,
                                pb_t_manafest_ps_trbs
                               )),
                       pb_t_manafest_vs_trbs)

all_vs_trbs <- setdiff(pb_t_manafest_vs_trbs, 
                       unique(c(sc_tissue_tet_ps_trbs, 
                                bulk_tissue_dex_ps_trbs, 
                                pb_t_manafest_ps_trbs
                               )))

length(all_ps_trbs)
length(all_vs_trbs)

# all TCR betas identified by multimer and MANAFEST assays
all_trbs <- c(all_ps_trbs, all_vs_trbs)
length(all_trbs)

length(all_ps_trbs)
length(unique(c(sc_tissue_tet_ps_trbs_novs, bulk_tissue_dex_ps_trbs_novs, 
                        pb_t_manafest_ps_trbs#, pb_m_manafest_ps_trbs
               )))

all_ps_trbs_df <- data.frame(pat_TRB_aa=all_trbs,
                             patient = gsub("(MSK.*)\\..*", "\\1", all_trbs))

all_ps_trbs_df$tissue_tet <- all_ps_trbs_df$pat_TRB_aa %in% sc_tissue_tet_ps_trbs
all_ps_trbs_df$tissue_dex <- all_ps_trbs_df$pat_TRB_aa %in% bulk_tissue_dex_ps_trbs
all_ps_trbs_df$pb_t_manafest_ps <- all_ps_trbs_df$pat_TRB_aa %in% pb_t_manafest_ps_trbs
all_ps_trbs_df$pb_t_manafest_vs <- all_ps_trbs_df$pat_TRB_aa %in% pb_t_manafest_vs_trbs

all_ps_trbs_df$nmethods <- rowSums(all_ps_trbs_df[,-c(1:2)])


all_ps_trbs_df$ps_category <- ifelse(all_ps_trbs_df$pb_t_manafest_vs, "Viral-specific", 
                                     ifelse(all_ps_trbs_df$nmethods > 1, "Peptide-specific (high confidence)",
                                            "Peptide-specific (low confidence)"))

table(subset(all_ps_trbs_df, patient=="MSK1263")$ps_category, useNA="ifany")

# ==== overlap of clones from different empirical methods ====
# helper function for venn diagram of TCR clone specificity between empirical methods
library(eulerr)

ps_group_pal <- c("pb_m_manafest"="#9fc5e8", "pb_t_manafest"="#fd745f",  "pb_t_manafest_v"="#9fcb92",
                  "sc_tissue_orig"="#fff2cc", "tissue_dex"="grey70", "tissue_tet"="white")

tcr_overlap_datasets <- function(cur_pat, with_sc_orig=FALSE, with_m_manafest=FALSE, 
                                 with_viral=FALSE, filter_vs=FALSE) {
    
    if (filter_vs) {
        tcr_list <- list(tissue_tet=sc_tissue_tet_ps_trbs_novs[grepl(cur_pat, sc_tissue_tet_ps_trbs_novs)],
                     tissue_dex=bulk_tissue_dex_ps_trbs_novs[grepl(cur_pat, bulk_tissue_dex_ps_trbs_novs)],
                     pb_t_manafest=pb_t_manafest_ps_trbs[grepl(cur_pat, pb_t_manafest_ps_trbs)]
                    )
    } else {
        tcr_list <- list(tissue_tet=sc_tissue_tet_ps_trbs[grepl(cur_pat, sc_tissue_tet_ps_trbs)],
                     tissue_dex=bulk_tissue_dex_ps_trbs[grepl(cur_pat, bulk_tissue_dex_ps_trbs)],
                     pb_t_manafest=pb_t_manafest_ps_trbs[grepl(cur_pat, pb_t_manafest_ps_trbs)]
                    )
    }
    
    if (with_m_manafest) {
        tcr_list$pb_m_manafest <- pb_m_manafest_ps_trbs[grepl(cur_pat, pb_m_manafest_ps_trbs)]

    }
    
    if (with_viral) {
        if (filter_vs) {
            tcr_list$pb_t_manafest_v <- all_vs_trbs[grepl(cur_pat, all_vs_trbs)]
        } else {
            tcr_list$pb_t_manafest_v <- pb_t_manafest_vs_trbs[grepl(cur_pat, pb_t_manafest_vs_trbs)]
        }
    }
    
    if (with_sc_orig) {
        tcr_list$sc_tissue_orig <- unique(cfi$pat_TRB_aa)[grepl(cur_pat, unique(cfi$pat_TRB_aa))]
    }
    
    tcr_overlap <- venn(tcr_list)

    plot(tcr_overlap, quantities = T, 
         fills = list(fill=ps_group_pal[names(tcr_list)], alpha=0.8),
         main = cur_pat)
}

# make venn diagrams of clone overlap
all_patients <- c("MSK1263", "MSK1302", "MSK1344")

venn_list <- list()
venn_list$ps_v_venns <- lapply(all_patients, tcr_overlap_datasets, 
                     with_sc_orig=FALSE, with_viral=TRUE
                    )

venn_list$ps_v_venns_novs <- lapply(all_patients, tcr_overlap_datasets, 
                     with_sc_orig=FALSE, with_viral=TRUE, filter_vs=TRUE
                    )

venn_list$ps_mm_venns <- lapply(all_patients, tcr_overlap_datasets, 
                      with_sc_orig=FALSE, with_viral=FALSE, with_m_manafest=TRUE
                     )

venn_list$ps_sc_mm_venns <- lapply(all_patients, tcr_overlap_datasets, 
                      with_sc_orig=TRUE, with_viral=FALSE, with_m_manafest=TRUE
                     )

venn_list$ps_sc_venns <- lapply(all_patients, tcr_overlap_datasets, 
                      with_sc_orig=TRUE, with_viral=FALSE, with_m_manafest=FALSE
                     )

venn_list$ps_sc_venns_novs <- lapply(all_patients, tcr_overlap_datasets, 
                      with_sc_orig=TRUE, with_viral=FALSE, with_m_manafest=FALSE, filter_vs=TRUE
                     )

venn_list$ps_sc_v_venns <- lapply(all_patients, tcr_overlap_datasets, 
                      with_sc_orig=TRUE, with_viral=TRUE, with_m_manafest=FALSE, filter_vs=FALSE
                     )

venn_list$ps_sc_v_venns_novs <- lapply(all_patients, tcr_overlap_datasets, 
                      with_sc_orig=TRUE, with_viral=TRUE, with_m_manafest=FALSE, filter_vs=TRUE
                     )

# Fig 3C: TCR tumor specificity overlap with total CD3+ dataset
pdf("plots/ps_beta_sc_venn_diagram.pdf", width=14, height=6)
do.call(grid.arrange, c(venn_list$ps_sc_venns_novs, nrow=1))
dev.off()

# Fig S10F: TCR specificity overlap
pdf("plots/ps_beta_v_novs_venn_diagram.pdf", width=14, height=6)
do.call(grid.arrange, c(venn_list$ps_v_venns_novs, nrow=1))
dev.off()

pdf("plots/ps_beta_sc_v_novs_venn_diagram.pdf", width=14, height=6)
do.call(grid.arrange, c(venn_list$ps_sc_v_venns_novs, nrow=1))
dev.off()

do.call(grid.arrange, c(lapply(venn_list, function(x) x[[1]]), nrow=1))
do.call(grid.arrange, c(lapply(venn_list, function(x) x[[2]]), nrow=1))

# high confidence neopeptide-specific clones
ps_trbs <- setdiff(intersect(sc_tissue_tet_ps_trbs, bulk_tissue_dex_ps_trbs), pb_t_manafest_vs_trbs)
length(ps_trbs)

save(all_ps_trbs_df, ps_trbs, all_ps_trbs, all_vs_trbs, pb_t_manafest_vs_trbs, 
     file = "r_objects/specific_trbs.RData")

# ==== map TCR specificity to tissue CD3+ dataset ====
# high confidence peptide-specific clones: TRBs found in tissue_tet, tissue_dex (+/- others)
# low confidence peptide-specific clones: other TRBs that were present in any of the 4 assays
# not peptide-specific: all other TRBs, undetected by any assay
cfi <- add_metadata_seurat(cfi, all_ps_trbs_df[,c("pat_TRB_aa", "ps_category")], 
                           join_by = "pat_TRB_aa", overwrite = T)

cfi$ps_category[is.na(cfi$ps_category)] <- "Unknown specificity"

cfi$ps_category <- factor(cfi$ps_category,
                          levels = c("Peptide-specific (high confidence)", 
                                     "Peptide-specific (low confidence)", 
                                     "Viral-specific", "Unknown specificity"))

table(unique(cfi@meta.data[,c("patient","pat_TRB_aa","ps_category")])[,c("patient","ps_category")])

ps_pal <- c(brewer.pal(n = 3, "Reds")[c(3,2)], "#9fcb92", "grey60")
names(ps_pal) <- levels(cfi$ps_category)

# Fig 3D: tumor reactivity score by TCR specificity
ggplot(subset(cfi[[]], 
              patient == "MSK1263" & cd_cluster_clone %in% c("CD8", "mixed-CD8")), 
       aes(x=ps_category, y=tumor_reactivity)) +
    geom_boxplot(aes(fill=ps_category, color=ps_category), alpha=0.6) +
    stat_compare_means(method = "wilcox.test", tip.length = 0.02, label = "p.signif",
                       comparisons = list(c("Peptide-specific (high confidence)", "Peptide-specific (low confidence)"),
                                          c("Peptide-specific (high confidence)", "Viral-specific"),
                                          c("Peptide-specific (low confidence)", "Viral-specific"),
                                          c("Peptide-specific (low confidence)", "Unknown specificity"),
                                          c("Peptide-specific (high confidence)", "Unknown specificity")
                                         ), 
                       ) +
    scale_fill_manual(values = ps_pal, guide='none') +
    scale_color_manual(values = ps_pal, guide='none') +
    labs(x="TCR specificity", y= "CD8 tumor reactivity score") +
    theme_bw() + expand_y_axis + rotatex + nobg + theme(plot.margin = unit(c(0,0,0,2.5), "cm"))

ggsave("plots/ps_category_tr_score_boxplot_withviral_psignif_msk1263.pdf", height=5, width=3)

clone_composition <- unique(data.frame(tmp %>%
    group_by(pat_TRB_aa, pat_cdr3, patient, ps_category) %>%
    mutate(clone_size = n()) %>%
    group_by(pat_TRB_aa, pat_cdr3, patient, ps_category, clone_size, pheno_cluster, .drop = F) %>%
    summarize(ncells = n()) %>%
    mutate(freq = ncells/clone_size) %>%
    arrange(desc(clone_size)) %>%
    subset(clone_size > 1 & grepl("CD8", pheno_cluster))))

clone_composition_region <- unique(data.frame(tmp %>%
    group_by(pat_TRB_aa, pat_cdr3, patient, ps_category) %>%
    mutate(clone_size = n()) %>%
    group_by(pat_TRB_aa, pat_cdr3, patient, ps_category, clone_size, tumor_status, .drop = F) %>%
    summarize(ncells = n()) %>%
    mutate(freq = ncells/clone_size) %>%
    arrange(desc(clone_size)) %>%
    subset(clone_size > 1)))

# Fig 3E: cluster frequency per TCR specificity
cur_pat <- "MSK1263"
min_clone_size <- 10

ggplot(subset(clone_composition,
              clone_size > min_clone_size &
              patient == cur_pat & ! pheno_cluster %in% c("CD8-NAIVE", "CD8-TCF1")),
       aes(x=ps_category, y=freq)) +
    geom_boxplot(aes(fill=ps_category, color=ps_category), alpha=0.6) +
    scale_fill_manual(values=ps_pal, limits=force, guide='none') +
    scale_color_manual(values=ps_pal, limits=force, guide='none') +
    facet_row(vars(pheno_cluster), space="free", scales="free") +
    stat_compare_means(method = "wilcox.test", label="p.signif",
                       comparisons = list(c("Peptide-specific (high confidence)", "Peptide-specific (low confidence)"),
                                          c("Peptide-specific (high confidence)", "Viral-specific"),
                                          c("Peptide-specific (low confidence)", "Viral-specific"),
                                          c("Peptide-specific (high confidence)", "Unknown specificity")                                         
                                         )
                      ) +
    labs(x="Peptide specificity", y="Frequency of cluster per clone", title=cur_pat) + 
    theme_bw() + 
    theme(legend.position="bottom", plot.margin=unit(c(0,0,0,3),"cm")) + nobg + rotatex + expand_y_axis

ggsave("plots/ps_category_pheno_freq_msk1263_psignif.pdf", height=5, width=10, useDingbats=F)

# helper function to make clonal bar plots
plot_clone_comps <- function(cur_ps_cat, clone_col, color_by, cur_pats, exclude_regions=NA) {

    # top clones expanded > 10 cells
    top_cur_trbs <- data.frame(top_expanded_clones %>% 
                             filter(patient %in% cur_pats & clone_size >= 10 & 
                                    ps_category == cur_ps_cat) %>%
                               top_n(wt = clone_size, n=70)
                              )
    top_cur_trbs[[clone_col]] <- factor(top_cur_trbs[[clone_col]], levels=top_cur_trbs[[clone_col]])

    # subset to cells in top clones
    tmp <- subset(cfi[[]], 
                  ! region_type %in% exclude_regions &
                  ps_category == cur_ps_cat & 
                  cd_cluster == "CD8" &
                  patient %in% cur_pats)
    tmp <- tmp[tmp[[clone_col]] %in% top_cur_trbs[[clone_col]], ]
    
    stopifnot(length(intersect(unique(tmp[[clone_col]]),
                               unique(top_cur_trbs[[clone_col]]))) == length(unique(tmp[[clone_col]])))
    tmp[[clone_col]] <- factor(tmp[[clone_col]], levels=top_cur_trbs[[clone_col]])

    # plot bar plot
    p <- ggplot() +
        geom_bar(data=tmp,
                 aes_string(x=clone_col, fill=color_by), position="fill") +
        geom_text(data=top_cur_trbs, aes_string(x=clone_col, label="clone_size"), y=1.04, size=2) +
        ggforce::facet_row(vars(ps_category)) +
        coord_cartesian(ylim=c(0,1.1), expand=F) +
        labs(x=paste0("TCR clone (", clone_col, ")"), y="# Cells", title=cur_pat) +
    guides(fill=guide_legend(nrow=1)) +
        theme_bw() + theme(axis.text.x=element_blank(), legend.position="bottom") + nobg
    
    if (color_by == "pheno_cluster") {
        p + scale_fill_manual(values=cluster_pal, limits=force, name="Cluster")
    } else if (color_by == "tumor_status") {
        p + scale_fill_manual(values=tumor_status_pal, limits=force, name="Region type")
    }
}

# Fig 3F-G: phenotype and region type of clones by TCR specificity
cur_pat <- "MSK1263"
cur_ps_cat <- "Peptide-specific (high confidence)"

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "pheno_cluster",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )
ggsave(file="plots/ps_clone_phenotypes_barplot_msk1263_highconf.pdf", width=12, height=4)

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "tumor_status",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )

ggsave(file="plots/ps_clone_tumorstatus_barplot_msk1263_highconf.pdf", width=12, height=4)

cur_ps_cat <- "Peptide-specific (low confidence)"

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "pheno_cluster",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )
ggsave(file="plots/ps_clone_phenotypes_barplot_msk1263_lowconf.pdf", width=12, height=4)

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "tumor_status",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )

ggsave(file="plots/ps_clone_tumorstatus_barplot_msk1263_lowconf.pdf", width=12, height=4)

cur_ps_cat <- "Unknown specificity"

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "pheno_cluster",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )
ggsave(file="plots/ps_clone_phenotypes_barplot_msk1263_notps.pdf", width=12, height=4)

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "tumor_status",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )

ggsave(file="plots/ps_clone_tumorstatus_barplot_msk1263_notps.pdf", width=12, height=4)

cur_ps_cat <- "Viral-specific"

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "pheno_cluster",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )
ggsave(file="plots/ps_clone_phenotypes_barplot_msk1263_viral.pdf", width=2, height=4)

plot_clone_comps(cur_ps_cat,
                 clone_col = "pat_cdr3",
                 color_by = "tumor_status",
                 cur_pats = cur_pat,
                 exclude_regions = c("Adrenal")
                )

ggsave(file="plots/ps_clone_tumorstatus_barplot_msk1263_viral.pdf", width=2, height=4)

# Fig 3J: TCR regionality by specificity
cur_df <- unique(subset(cfi[[]], 
                        ! is.na(regionality_tumor) & clone_size > 1 &
                        patient != "MSK1344" &  cd_clone %in% c("CD8", "mixed-CD8")
                       )[,c("regionality", "regionality_tumor", "pat_TRB_aa", "pat_cdr3", 
                            "patient", "ps_category")])
clone_counts_df <- cur_df %>% group_by(ps_category) %>% summarize(nclones=n())

ggplot(cur_df, 
       aes(x=ps_category)) +
    geom_bar(aes(fill=regionality_tumor), position="fill") +
    geom_text(data=clone_counts_df, aes(label=nclones), size=2, y=1.04) +
    scale_fill_manual(values = regionality_tumor_pal) +
    coord_cartesian(ylim=c(0,1.1), expand=F) +
    labs(x="Peptide specificity") +
    theme_bw() + rotatex + nobg

ggsave("plots/ps_category_tcr_regionality_tumor_bar_expanded.pdf", height=4.5, width=4)

# concordance of empirical vs. transcriptional tumor reactivity
cd8_tumor_clone_df <- readRDS(file="r_objects/cd8_tumor_clone_df.rds")

head(cd8_tumor_clone_df)

ps_tr_concordance <- dplyr::left_join(subset(cfi@meta.data, patient=="MSK1263"),
                                      cd8_tumor_clone_df[, c("pat_cdr3", "tr_status")],
                                      by = "pat_cdr3"
                                     ) %>%
    filter(!is.na(tr_status)) %>%
    group_by(ps_category) %>% mutate(ncells_ps=n()) %>%
    group_by(tr_status, ps_category, ncells_ps) %>% summarize(ncells=n()) %>%
    mutate(prop=ncells/ncells_ps)
                                      
# Fig S10G: concordance heatmap
ggplot(subset(ps_tr_concordance, !is.na(ps_tr_concordance$tr_status)),
       aes(x=tr_status, y=ps_category)) +
    geom_tile(aes(fill=prop)) +
    geom_text(aes(label=ncells), color='white') +
    scale_fill_viridis(option = "viridis", name="Proportion of cells") +
    coord_cartesian(expand=F) +
    labs(x="Tumor reactivity", y="Peptide specificity") +
    theme_bw() + nobg

ggsave("plots/ps_tr_concordance_msk1263.pdf", width=5.5, height=3)

# ==== regional DE of neopeptide-specific clones ====
ps_hc_cells <- subset(cfi, 
              ps_category == "Peptide-specific (high confidence)" &
              patient != "MSK1344" & cd_cluster_clone %in% c("CD8", "mixed-CD8"))

table(ps_hc_cells[[c("patient", "tumor_status")]])

# LN vs. tumor DE
DefaultAssay(ps_hc_cells) <- "RNA"

Idents(ps_hc_cells) <- "region_type"
ident1 = "Tumor"; ident2 = "LN"

h_tcr_genes <- get_tcr_genes(ps_hc_cells, "human")

cur_de <- FindMarkers(subset(ps_hc_cells, patient=="MSK1263"), max.cells.per.ident = 500, 
                           assay = "RNA", `ident.1` = ident1, `ident.2`= ident2,
                           features = setdiff(rownames(ps_hc_cells), h_tcr_genes))

cur_de$gene <- rownames(cur_de)
cur_de$cluster <- ifelse(cur_de$avg_log2FC > 0, ident1, ident2)

cur_de <- data.frame(cur_de %>% 
                         arrange(desc(avg_log2FC), p_val_adj) %>% 
                         group_by(cluster) %>% 
                         mutate(fc_rank = rank(dplyr::desc(abs(avg_log2FC)), ties.method = "first"),
                                p_rank = rank(p_val_adj, ties.method = "first")
                               ))

# Fig S10H: tumor vs. LN DE
make_de_plot(cur_de, paste(ident1, "vs.", ident2),
             "Peptide-specific (high-confidence) clones: MSK1263", "avg_log2FC", 
             fc_rank_cutoff = 20, p_rank_cutoff = 20) + 
    scale_color_manual(values=region_type_pal, name = "Region type", limits=force)  +
    coord_cartesian(ylim=c(0,20), xlim=c(-2.5,2.5)) 

ggsave("plots/ps_ln_tumor_de_volcano.pdf", height=7, width=7, useDingbats=F)

# VT vs. NVT DE
DefaultAssay(ps_hc_cells) <- "RNA"

Idents(ps_hc_cells) <- "tumor_status"
ident1 = "Viable Tumor"; ident2 = "No Viable Tumor"

h_tcr_genes <- get_tcr_genes(ps_hc_cells, "human")

cur_de <- FindMarkers(subset(ps_hc_cells, patient=="MSK1263"),
                           assay = "RNA", `ident.1` = ident1, `ident.2`= ident2,
                           features = setdiff(rownames(ps_hc_cells), h_tcr_genes))

cur_de$gene <- rownames(cur_de)
cur_de$cluster <- ifelse(cur_de$avg_log2FC > 0, ident1, ident2)

cur_de <- data.frame(cur_de %>% 
                         arrange(desc(avg_log2FC), p_val_adj) %>% 
                         group_by(cluster) %>% 
                         mutate(fc_rank = rank(dplyr::desc(abs(avg_log2FC)), ties.method = "first"),
                                p_rank = rank(p_val_adj, ties.method = "first")
                               ))

# Fig S10I: viable tumor vs. no viable tumor DE
make_de_plot(cur_de, paste(ident1, "vs.", ident2),
             "Peptide-specific (high-confidence) clones: MSK1263", "avg_log2FC", 
             fc_rank_cutoff = 20, p_rank_cutoff = 20) + 
    scale_color_manual(values=tumor_status_pal, name = "Region type", limits=force)  +
    coord_cartesian(ylim=c(0,60), xlim=c(-1,1))

ggsave("plots/ps_vt_nvt_de_volcano.pdf", height=7, width=7, useDingbats=F)

# ==== peripheral persistence ====
load(file="r_objects/bulk_data.RData")
load(file="r_objects/bulk_revision_data_obj.RData")

head(unique(cfi@meta.data[,c("pat_TRB_aa", "high_conf_ps_trb")]))

bulk_tcr_bytimepoint <- data.frame(bulk_tcr_bytimepoint %>% 
                                   group_by(pat_TRB_aa) %>% 
                                   mutate(beta_clone_size=sum(sum_templates),
                                          beta_clone_freq=sum(sum_freq)))

bulk_tcr_bytimepoint$scaled_size <- bulk_tcr_bytimepoint$sum_templates / bulk_tcr_bytimepoint$beta_clone_size
bulk_tcr_bytimepoint$scaled_freq <- bulk_tcr_bytimepoint$sum_freq / bulk_tcr_bytimepoint$beta_clone_freq

ntimepoints_per_beta <- data.frame(bulk_tcr_bytimepoint %>% group_by(patient, pat_TRB_aa) %>% 
                                   summarize(ntimepoints=n_distinct(date)))

max_timepoints <- data.frame(bulk_tcr_bytimepoint %>% group_by(patient) %>%
                             summarize(total_timepoints=n_distinct(date)))
rownames(max_timepoints) <- max_timepoints$patient
max_timepoints$patient <- NULL
max_timepoints

bulk_tcr_bytimepoint <- dplyr::left_join(bulk_tcr_bytimepoint,
                unique(cfi@meta.data[,c("pat_TRB_aa", "ps_category")]), by="pat_TRB_aa")

table(subset(cfi@meta.data, cd_cluster == "CD8")$ps_category, useNA="ifany")
prop.table(table(subset(cfi@meta.data, cd_cluster == "CD8")$ps_category, useNA="ifany"))

# Fig 4D: longitudinal persistence of clones
pat <- "MSK1263"

all_tp_betas <- subset(ntimepoints_per_beta, patient==pat & 
                       ntimepoints == max_timepoints[pat,"total_timepoints"])$pat_TRB_aa

ggplot(subset(bulk_tcr_bytimepoint, patient==pat & pat_TRB_aa %in% all_tp_betas), 
       aes(x=date, y=scaled_size, color=ps_category)) +
    rasterize(geom_point(), dpi=600) + 
    rasterize(geom_line(aes(group=pat_TRB_aa), alpha=0.3), dpi=600) +
    facet_grid(ps_category~patient) +
    stat_summary(fun="mean", aes(group=1), geom="line", colour="black") +
    stat_summary(geom="ribbon", 
                 fun.data = mean_cl_normal, fun.args=list(conf.int=0.95),
                 fill="grey", color=NA, alpha=0.3) +#     scale_y_log10() +
    scale_color_manual(values=ps_pal, guide='none') +
    coord_cartesian(ylim=c(0,1)) +
    labs(x="Timepoint", y="Bulk TCRb frequency (scaled per clone)") +
    scale_x_date(date_labels = "%m/%d/%Y", 
                 breaks = unique(subset(bulk_tcr_bytimepoint, patient==pat)$date)
                ) + 
    theme_bw() + rotatex 

ggsave("plots/ps_category_persistence_msk1263_rasterized.pdf", width=2.5, height=8.2)

save(bulk_tcr_bytimepoint, ntimepoints_per_beta, max_timepoints, 
    file="r_objects/bulk_revision_data_obj.RData")