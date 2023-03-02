# ==== set up ====
packages <- c("Seurat", "magrittr", "dplyr", "reshape2", "pheatmap",
              "ggplot2", "ggrepel", "ggsci", "ggpubr", "ggrastr",
              "RColorBrewer", "wesanderson", "viridis",
              "gridExtra", "grid", "cowplot")

invisible(lapply(packages, library, character.only = TRUE))

source("~/oak/scripts/sc_utils.R")

# load data
cfi <- readRDS(file="r_objects/cfi.rds")
load(file="r_objects/all_color_pals.RData")
load(file="r_objects/regionality_cd.RData") # tcr_regionality_df

# ==== classify regions as tumor, normal, or LN ====
region_classification <- c(rep("Adrenal", 3), "LN", "Normal", rep("Tumor", 8),
                           "LN", "Normal", rep("Tumor", 8),
                           rep("Tumor", 7), "Normal")
names(region_classification) <- sort(unique(cfi$orig.ident))

# pathology data
pathology_data <- read.csv("region_pathology_data.csv", stringsAsFactors = F)
rownames(pathology_data) <- pathology_data$region

add_region_metadata <- function(df=NULL, sample_column=NULL, region_type=TRUE, patient=TRUE, path=TRUE) {
    df[,sample_column] <- as.character(df[,sample_column])
    
    if(region_type) { 
        # add region type
        df$region_type <- region_classification[df[,sample_column]]   
    }
    if(patient) {
        # add patient information
        df$patient <- gsub("^(MSK\\w+)_.*", "\\1", df[,sample_column])
    }
    if(path) {
        df$tumor_status <- pathology_data[df[,sample_column], "tumor_status"]
        df$percent_viable_tumor <- pathology_data[df[,sample_column], "percent_viable_tumor"]
    }
    return(df)
}

pathology_data <- add_region_metadata(pathology_data, "region", region_type=T, patient=F)
cfi@meta.data <- add_region_metadata(cfi@meta.data, "orig.ident", region_type=F, patient=F, path=TRUE)

# ==== CD3+ data overview ====
# Fig 1C: phenotype cluster UMAP
DimPlot(cfi, group.by = "pheno_cluster") + scale_color_manual(values=cluster_pal)

# Fig S3B: phenotype cluster UMAP by sample
DimPlot(cfi, group.by = "pheno_cluster", split.by="orig.ident", ncol=8) + 
    scale_color_manual(values=cluster_pal) +
    theme_classic()

# clonal expansion
cfi$pat_cdr3 <- paste(cfi$patient, cfi$cdr3s_nt, sep=".") # clone ID

cfi_tcrab_cells <- subset(cfi[[]], chains_captured=="AB")

cfi_tcrab_clones <- data.frame(cfi_tcrab_cells %>% group_by(pat_cdr3) %>% summarize(clone_size=n()))

tmp <- dplyr::left_join(cfi@meta.data, cfi_tcrab_clones, by="pat_cdr3")
stopifnot(identical(tmp$pat_cdr3, cfi@meta.data$pat_cdr3))
rownames(tmp) <- rownames(cfi[[]])
cfi@meta.data <- tmp

cfi$clone_size_bin <- ifelse(cfi$clone_size == 1, "singleton", 
                             ifelse(cfi$clone_size <= 10, "1-10", 
                                    ifelse(cfi$clone_size <= 50, "11-50", 
                                           ifelse(cfi$clone_size <= 100, "51-100", 
                                                  ifelse(cfi$clone_size > 100, ">100", "no tcr")
                                          ))
                             )
                            )

cfi$clone_size_bin <- factor(cfi$clone_size_bin,
                             levels = c("singleton", "1-10", "11-50", "51-100", ">100"))

clone_size_pal <- viridis(length(levels(cfi$clone_size_bin)))

# Fig 1E: clone size UMAP
DimPlot(subset(cfi, clone_size_bin != "no tcr"), pt.size = 3,
        group.by = "clone_size_bin", cols = clone_size_pal, raster=T)

ggsave("plots/Fig1E_clone_size_bin_umap.pdf", height=5, width=6.2)

FeaturePlot(cfi_tcrab_cells, cols = c("grey", "red"), reduction="umap", 
            features = "clone_size", min.cutoff = "q05", max.cutoff="q95") +
    scale_color_viridis(trans="log10") +
    labs(title="Clone sizes")

# cluster markers
Idents(cfi) <- "pheno_cluster"
cluster.averages <- AverageExpression(cfi, return.seurat=TRUE) 

top_genes <- cluster_markers %>% arrange(cluster_annotation) %>%
    group_by(cluster_annotation) %>% top_n(n = 10, wt=avg_logFC)

# select T cell markers
features <- c("PDCD1", "CTLA4", "TIGIT", "HAVCR2", "LAG3", "ENTPD1", "TOX",
              "FOXP3", "TNFRSF9", "TNFRSF18", "TNFRSF4",
              "CCR7", "SELL", "IL7R", "LEF1", "TCF7",
              "PRDM1", "CXCR6", "ICOS", "CXCL13", 
              "BCL6", "CXCR5", "CD4", "CD8A",
              "GZMK", "CCL4", "GZMH",
              "PRF1", "GZMA", "GZMB", "NKG7", "IFNG", "KLRG1", "KLRD1", 
              "CD69", "ITGAE", "ZNF683", "CXCR3",
              "MKI67", "TUBB", "CD38", "CD274"
             )

# Fig 1D: cluster marker heatmap
pheatmap(cluster.averages@assays$RNA[features,],
         color = colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100),
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         scale = "row",
         angle_col=45,
         filename="plots/cluster_markers_heatmap.pdf", height=8.5, width=5
)

# cluster proportions
cluster_proportions_bar <- function(seurat_obj, cluster_var, group_var) {
    cluster_counts_by_group <- sapply(as.character(sort(unique(seurat_obj@meta.data[,group_var]))), 
                                      function(x) { 
                                          cur_cells <- seurat_obj@meta.data[seurat_obj@meta.data[,group_var]==x,]
                                          table(cur_cells[,cluster_var])}
                                     )
    cluster_props_by_group <- prop.table(cluster_counts_by_group, 2)

    cluster_props_df <- melt(cluster_props_by_group, direction = "long",
                         varnames = c(cluster_var, group_var), value.name = "proportion")
    
    group_levels <- levels(cluster_props_df[,group_var])
    cluster_props_df <- add_region_metadata(cluster_props_df, group_var)
    cluster_props_df[,group_var] <- factor(cluster_props_df[,group_var], levels=group_levels)
    
    # generate barplot
    p <- ggplot(cluster_props_df, aes_string(x=group_var, y="proportion", fill=cluster_var)) +
        geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + coord_cartesian(expand=F) +
        labs(y="proportion of cells") + theme_bw() +
        theme(axis.text.x=element_text(angle=45, hjust=1))
    
    clust_prop_res <- list(proportions=cluster_props_df, plot=p)
    return(clust_prop_res)
}

clust_props <- cluster_proportions_bar(cfi, "pheno_cluster", "orig.ident")

clust_props$proportions$cd_cluster <- gsub("^(CD[48])-.*", "\\1", clust_props$proportions$pheno_cluster)
clust_props$proportions$orig.ident <- as.character(clust_props$proportions$orig.ident)
clust_props$proportions <- add_region_metadata(clust_props$proportions, "orig.ident", path=T)
clust_props$proportions$orig.ident <- factor(clust_props$proportions$orig.ident, levels=sample_order)

# Fig S3A: CD8 cluster representation per sample
ggplot(subset(clust_props$proportions, cd_cluster=="CD8"), 
       aes(x=orig.ident, y=proportion, fill=pheno_cluster)) +
    geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + 
    coord_cartesian(expand=F) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    scale_fill_manual(values=cluster_pal, name="Cluster", drop=T, limits = force) + 
    facet_grid(.~patient, scales="free", space="free") +
    labs(x="Region", y="Proportion of cells in CD8 clusters")

ggsave("plots/clust_prop_cd8.pdf", height=3.5, width=11, useDingbats=F)

# Fig S3B: CD4 cluster representation per sample
ggplot(subset(clust_props$proportions, cd_cluster=="CD4"), 
       aes(x=orig.ident, y=proportion, fill=pheno_cluster)) +
    geom_bar(stat="identity", position = position_fill(reverse = TRUE)) + 
    coord_cartesian(expand=F) +
    theme_bw() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) + 
    scale_fill_manual(values=cluster_pal, name="Cluster", drop=T, limits = force) + 
    facet_grid(.~patient, scales="free", space="free") +
    labs(x="Region", y="Proportion of cells in CD4 clusters")

ggsave("plots/clust_prop_cd4.pdf", height=3.5, width=11, useDingbats=F)

pheno_props <- subset(cfi[[]], 
                      patient != "MSK1344" & region_type != "Adrenal") %>%
    group_by(tumor_status) %>% mutate(ncells_total=n()) %>%
    group_by(tumor_status, ncells_total, pheno_cluster) %>% 
    summarize(ncells=n()) %>%
    mutate(prop=ncells/ncells_total)

pheno_prop_wide <- dcast(pheno_props[,c("tumor_status", "pheno_cluster", "prop")],
                         tumor_status~pheno_cluster)
rownames(pheno_prop_wide) <- pheno_prop_wide$tumor_status
pheno_prop_wide$tumor_status <- NULL

cd4_clusters <- colnames(pheno_prop_wide)[grepl("CD4", colnames(pheno_prop_wide))]
cd8_clusters <- colnames(pheno_prop_wide)[grepl("CD8", colnames(pheno_prop_wide))]

# Fig 1F: CD8 cluster frequency per region type
pheatmap::pheatmap(t(pheno_prop_wide[,cd8_clusters]), scale="row", color = magma(100), 
                   cluster_cols = FALSE,
                   treeheight_row = 15, treeheight_col = 15, angle=45, 
                   filename="plots/cluster_heatmap_cd8.pdf", width=3, height=2.5,
                  )

# Fig 1G: CD4 cluster frequency per region type
pheatmap::pheatmap(t(pheno_prop_wide[,cd4_clusters]), scale="row", color = magma(100), 
                   cluster_cols = FALSE,
                   treeheight_row = 15, treeheight_col = 15, angle=45, 
                   filename="plots/cluster_heatmap_cd4.pdf", width=2.7, height=2.5,
                  )

# ==== classify TCR clones as CD4 or CD8 based on CD4/8 cluster membership ====
cfi$cd_cluster <- gsub("^(CD[48])-.*", "\\1", cfi$pheno_cluster)
complete_tcr_cells <- subset(cfi[[]], chains_captured == "AB")

clone_cd_cluster <- data.frame(complete_tcr_cells %>% 
                               group_by(pat_cdr3, cd_cluster) %>% 
                               summarize(ncells=n()))

clone_cd_cluster_wide <- dcast(clone_cd_cluster, pat_cdr3 ~ cd_cluster, value.var="ncells")
clone_cd_cluster_wide[is.na(clone_cd_cluster_wide)] <- 0

clone_cd_cluster_wide$clone_size <- rowSums(clone_cd_cluster_wide[,c(2:4)])
clone_cd_cluster_wide$CD4_percent <- clone_cd_cluster_wide$CD4 / clone_cd_cluster_wide$clone_size
clone_cd_cluster_wide$CD8_percent <- clone_cd_cluster_wide$CD8 / clone_cd_cluster_wide$clone_size

clone_cd_cluster_wide$cd_cluster_clone <- ifelse(clone_cd_cluster_wide$CD4_percent == 1, "CD4", 
                               ifelse(clone_cd_cluster_wide$CD8_percent == 1, "CD8",
                                     ifelse(clone_cd_cluster_wide$CD8_percent >= 0.75, "mixed-CD8", 
                                            ifelse(clone_cd_cluster_wide$CD4_percent >= 0.75, "mixed-CD4",
                                                   ifelse(clone_cd_cluster_wide$CD8 + clone_cd_cluster_wide$CD4 == 0, "MAIT",
                                                          "mixed")))))

clone_cd_cluster_wide$cd_cluster_clone <- factor(clone_cd_cluster_wide$cd_cluster_clone, 
                                                 levels=c("CD4", "CD8", "MAIT", "mixed-CD4", "mixed-CD8", "mixed")
                                                )

cd_cluster_clone_counts <- table(clone_cd_cluster_wide$cd_cluster_clone)

clone_cd_cluster_wide$cd_cluster_clone_facet <- paste0(clone_cd_cluster_wide$cd_cluster_clone,
                                                       " (", 
                                                       cd_cluster_clone_counts[clone_cd_cluster_wide$cd_cluster_clone], 
                                                       " clones)")

clone_cd_cluster_wide$cd_cluster_clone_facet <- factor(clone_cd_cluster_wide$cd_cluster_clone_facet,
                                                       levels=sort(unique(clone_cd_cluster_wide$cd_cluster_clone_facet))[c(1,2,3,5,6,4)]
                                                      )

prop.table(cd_cluster_clone_counts)*100

clone_cd_cluster_wide$cd_category <- ifelse(clone_cd_cluster_wide$cd_cluster_clone %in% c("CD4", "mixed-CD4"), 
                                            "CD4",
                                            ifelse(clone_cd_cluster_wide$cd_cluster_clone %in% c("CD8", "mixed-CD8"),
                                                   "CD8",
                                                   ifelse(clone_cd_cluster_wide$cd_cluster_clone == "MAIT", 
                                                          "MAIT", 
                                                          "mixed")))

clone_cd_cluster_wide$cd_category <- factor(clone_cd_cluster_wide$cd_category, 
                                            levels=c("CD4", "CD8", "mixed", "MAIT"))

ggplot(clone_cd_cluster_wide, aes(x=CD8_percent, y=clone_size)) +
    geom_point(alpha=0.5, aes(color=CD4_percent), size=2) +
    facet_wrap(~cd_cluster_clone_facet) +
    scale_color_viridis(name="% CD4 cells in clone") +
    scale_y_log10() +
    labs(x="% CD8 cells in clone", y="Clone size") +
    theme_bw()

ggsave("plots/clone_cd_scatter_cd8pct_clonesize.pdf", width=8, height=4)

ggplot(clone_cd_cluster_wide, aes(x=CD4_percent, y=clone_size)) +
    geom_point(alpha=0.5, aes(color=CD8_percent), size=2) +
    facet_wrap(~cd_cluster_clone_facet) +
    scale_color_viridis(name="% CD8 cells in clone") +
    scale_y_log10() +
    labs(x="% CD4 cells in clone", y="Clone size") +
    theme_bw()

ggsave("plots/clone_cd_scatter_cd4pct_clonesize.pdf", width=8, height=4)

# Fig S4B: clone CD categorization
p <- ggplot(clone_cd_cluster_wide, aes(x=CD4_percent, y=clone_size)) +
    geom_point(alpha=0.5, aes(color=CD8_percent), size=2) +
    facet_wrap(~cd_category, nrow=2) +
    scale_color_viridis(name="% CD8 cells in clone") +
    scale_y_log10() +
    labs(x="% CD4 cells in clone", y="Clone size") +
    theme_bw() + theme(legend.position="top")

rasterize(p, layers='Point', dpi=500)
ggsave("plots/clone_cd_cat_scatter_cd8pct_clonesize.pdf", width=6, height=4)

tmp <- dplyr::left_join(cfi@meta.data, 
                        clone_cd_cluster_wide[,c("pat_cdr3", "cd_cluster_clone")], 
                        by = "pat_cdr3")
stopifnot(identical(tmp$sample_cell_id, cfi@meta.data$sample_cell_id))
rownames(tmp) <- rownames(cfi@meta.data)
cfi@meta.data <- tmp

cfi$cd_category <- ifelse(cfi$cd_cluster_clone %in% c("CD4", "mixed-CD4"), "CD4",
                          ifelse(cfi$cd_cluster_clone %in% c("CD8", "mixed-CD8"), "CD8",
                                 ifelse(cfi$cd_cluster_clone == "MAIT", "MAIT", "mixed")))

cfi$cd_category <- factor(cfi$cd_category, levels= c("CD4", "CD8", "mixed", "MAIT"))

counts_df <- subset(cfi[[]], !is.na(cd_cluster_clone)) %>% 
    group_by(cd_cluster_clone, cd_category) %>% summarize(ncells=n())

# Fig S4C: phenotypes of cells within each clone CD category
ggplot(subset(cfi[[]], !is.na(cd_cluster_clone)), aes(x=cd_cluster_clone)) +
    geom_bar(aes(fill=pheno_cluster), position="stack") +
    geom_text(data=counts_df, aes(label=ncells, y=ncells+1500)) +
    scale_fill_manual(values=cluster_pal, name="Cluster") + 
    facet_grid(~cd_category, scales="free_x", space="free_x") +
    coord_cartesian(expand=F, ylim=c(0, 70000)) +
    labs(x="Clone CD type", y="# Cells") +
    theme_bw() + rotatex

ggsave("plots/clone_cd_cluster_barplot_stack.pdf", width=6, height=5)

ggplot(subset(cfi[[]], !is.na(cd_cluster_clone)), aes(x=cd_cluster_clone)) +
    geom_bar(aes(fill=pheno_cluster), position="fill") +
    geom_text(data=counts_df, aes(label=ncells), y=1.02) +
    scale_fill_manual(values=cluster_pal, name="Cluster") +
    facet_grid(~cd_category, scales="free_x", space="free_x") +
    coord_cartesian(expand=F, ylim=c(0, 1.05)) +
    labs(x="Clone CD type", y="Proportion of cells") +
    theme_bw() + rotatex

ggsave("plots/clone_cd_cluster_barplot_fill.pdf", width=6, height=5)

# ==== gene module scoring ====
# transcriptional signatures from literature (Table S2)
genesets <- list()
genesets$dysfunction <- scan("genesets/li_dysfunctional_geneset.txt", what = "character")
genesets$tumor_reactivity <- c("ENTPD1", "ITGAE", "PDCD1", "TNFRSF9", "TNFRSF4", "CXCL13", "MKI67",
                               "HAVCR2", "LAG3")

guegen_genes <- read.table("genesets/Guegen_abd5778_Table_S3.csv", 
                           sep=",", stringsAsFactors=F, header=T)
colnames(guegen_genes) <- c("exhausted_progenitor", "terminal_exhausted")
progen_genes_cleaned <- guegen_genes$exhausted_progenitor[!grepl("^MT-", guegen_genes$exhausted_progenitor) &
                                                          !grepl("^RP[SL]", guegen_genes$exhausted_progenitor)]
genesets$progenitor_cleaned <- progen_genes_cleaned

genesets$oliviera_virus_specific <- scan(file="genesets/oliviera_wu_virus_specific.csv",
                                        what = "character")
genesets$caushi_influenza_specific <- scan(file="genesets/caushi_smith_influenza_specific.csv",
                                        what = "character")
genesets$caushi_MANA_specific <- scan(file="genesets/caushi_smith_MANA_specific.csv",
                                        what = "character")
genesets$oliviera_tumor_specific <- scan(file="genesets/oliviera_wu_tumor_specific.csv",
                                        what = "character")

genesets$neotcr_cd4_40 <- scan("/oak/stanford/groups/satpathy/users/joypai/genesets/neotcr_cd4_40.csv", 
                              what = "character", skip = 1)
genesets$neotcr_cd8_all <- scan("/oak/stanford/groups/satpathy/users/joypai/genesets/neotcr_cd8_all.csv", 
                              what = "character", skip = 1)

saveRDS(genesets, file="r_objects/genesets.rds")

# add gene module scores per cell
DefaultAssay(cfi) <- "RNA"

cfi <- AddModuleScore(cfi, features = genesets, name = names(genesets), assay="RNA")
cfi <- rename_module_cols(cfi, genesets)

# ==== tumor reactivity ====
# CD8 tumor reactivity
tr_pheno_order <- data.frame(cfi[[]] %>% group_by(pheno_cluster) %>% 
                             summarize(mean_tr=mean(tumor_reactivity)) %>% arrange(desc(mean_tr)))

# Fig S8C: CD8 tumor reactivity per cluster
ggplot(subset(cfi[[]], cd_cluster=="CD8"),
       aes(x=factor(pheno_cluster, levels=tr_pheno_order$pheno_cluster), y=tumor_reactivity)) +
    geom_boxplot(aes(fill=pheno_cluster), outlier.shape=NA) +
    stat_compare_means(method = "wilcox.test", paired = F, label = "p.signif", 
                       step.increase = 0.1, tip.length = 0.01,
                       comparisons = list(c("CD8-EXH", "CD8-TRM"),
                                          c("CD8-EXH", "CD8-EFF"),
                                          c("CD8-EXH", "CD8-GZMK"),
                                          c("CD8-PROLIF-EXH", "CD8-EXH"),
                                          c("CD8-PROLIF-EXH", "CD8-TRM"),
                                          c("CD8-PROLIF-EXH", "CD8-EFF"),
                                          c("CD8-PROLIF-EXH", "CD8-GZMK")
                                         )) +
    scale_fill_manual(values=cluster_pal, guide=F) +
    labs(x="Cluster", y="Tumor reactivity score") +
    theme_classic() + rotatex

ggsave("plots/tr_score_pheno.pdf", height=5, width=3, useDingbats=F)

ggplot(subset(cfi[[]], cd_cluster=="CD8"),
       aes(x=factor(pheno_cluster, levels=tr_pheno_order$pheno_cluster), y=tumor_reactivity)) +
    geom_boxplot(aes(fill=pheno_cluster), outlier.shape=NA) +
    coord_cartesian(ylim=c(-0.5, 2)) +
    facet_grid(~patient, scales="free", space="free") +
    scale_fill_manual(values=cluster_pal, guide=F) +
    labs(x="Cluster", y="Tumor reactivity score") +
    theme_classic() + rotatex

# Fig S8D: CD8 tumor reactivity per region
ggplot(subset(cfi[[]], cd_cluster=="CD8"),
       aes(x=factor(tumor_status, levels=c("LN", "Normal", "No Viable Tumor", "Viable Tumor")), 
           y=tumor_reactivity)) +
    geom_boxplot(aes(fill=tumor_status), outlier.shape=NA) +
    coord_cartesian(ylim=c(-0.5, 2)) +
    stat_compare_means(label.y = c(1.7, 1.4, 1.1), 
                        comparisons = list(c("Viable Tumor", "LN"), 
                                          c("Viable Tumor", "Normal"),
                                          c("Viable Tumor", "No Viable Tumor")),
                       label="p.signif", tip.length=0.01
                      ) +
    scale_fill_manual(values=tumor_status_pal, guide=F) +
    labs(x="Region type", y="Tumor reactivity score") +
    theme_classic() + rotatex

ggsave("plots/tr_score_region.pdf", height=4, width=2.5, useDingbats=F)

# CD4 tumor reactivity
neotcr_cd4_pheno_order <- data.frame(cfi[[]] %>%
                                     group_by(pheno_cluster) %>% 
                                     summarize(mean_neotcr_cd4=mean(neotcr_cd4_40)) %>% 
                                     arrange(desc(mean_neotcr_cd4)))

# Fig S8F: CD4 tumor reactivity per cluster
ggplot(subset(cfi[[]], cd_cluster=="CD4"),
       aes(x=factor(pheno_cluster, levels=neotcr_cd4_pheno_order$pheno_cluster), y=neotcr_cd4_40)) +
    geom_boxplot(aes(fill=pheno_cluster), outlier.shape=NA) +
    stat_compare_means(method = "wilcox.test", paired = F, label = "p.signif", 
                       step.increase = 0.1, tip.length = 0.01,
                       comparisons = list(c("CD4-TFH2", "CD4-TREG"),
                                          c("CD4-TREG", "CD4-TFH1"),
                                          c("CD4-TREG", "CD4-EFF2"),
                                          c("CD4-TREG", "CD4-EFF1"),
                                          c("CD4-TFH2", "CD4-EFF2"),
                                          c("CD4-TFH2", "CD4-EFF1")
                                         )) +
    scale_fill_manual(values=cluster_pal, guide=F) +
    labs(x="Cluster", y="NeoTCR-CD4 score") +
    theme_classic() + rotatex

ggsave("plots/neotcr_cd4_score_pheno.pdf", height=4.2, width=2.8, useDingbats=F)

# classify cells as tumor reactive high or low
tr_boundary = 0

ggplot(cfi[[]], aes(x=tumor_reactivity)) +
    geom_histogram(binwidth=0.01) + 
    geom_vline(xintercept = tr_boundary, lty=2, color='red') +
    labs(x="Tumor reactivity score (per cell)", y="# of cells") +
    theme_bw()

cfi$cell_tr_category <- ifelse(cfi$tumor_reactivity >= tr_boundary, 
                               "TRhi", "TRlo")

cfi$cell_neotcr_cd8_category <- ifelse(cfi$neotcr_cd8_all >= neotcr_cd8_boundary,
                                       "neoTCR_CD8_hi", "neoTCR_CD8_lo")

# evaluate tumor reactivity signatures
tumor_specific_score_df <- cfi[[c("sample_cell_id", "tumor_reactivity", "neotcr_cd8_all",
                                  "cell_tr_category", "cell_neotcr_cd8_category", 
                                  "oliviera_tumor_specific", "caushi_MANA_specific", 
                                  "oliviera_virus_specific", "caushi_influenza_specific",
                                  "neotcr_cd8_all", "neotcr_cd4_40")]]

tumor_specific_score_long <- melt(tumor_specific_score_df, 
                                  id.vars = c("sample_cell_id", "tumor_reactivity", "cell_tr_category",
                                             "neotcr_cd8_all", "cell_neotcr_cd8_category"),
                                  variable.name = "gene_signature", value.name="score"
                                 )

head(tumor_specific_score_long)

# correlation of signatures
tumor_specific_sig_cor <- cor(tumor_specific_score_df[,c("tumor_reactivity", 
                                                         "oliviera_tumor_specific", 
                                                         "caushi_MANA_specific", 
                                                         "oliviera_virus_specific", 
                                                         "caushi_influenza_specific",
                                                         "neotcr_cd8_all")],
                              method = "pearson"
                             )

# Fig S8A: correlation of signatures
pheatmap(tumor_specific_sig_cor, angle=90, treeheight_row = 15, treeheight_col = 15,
         filename = "plots/tr_signature_corr.pdf", height=4.5, width=5
        )

# Fig S8B: tumor reactivity scores
gtr_pal <- brewer.pal(3, "Reds")[c(1,3)]
names(tr_pal) <- c("TRlo", "TRhi")

gplot(tumor_specific_score_long, 
       aes(x=cell_tr_category, y=score)) +
    geom_boxplot(aes(fill=cell_tr_category), outlier.shape=NA) +
    stat_compare_means(method = "t.test", label="p.signif", vjust = 0.5,
                       tip.length = 0.01,
                       comparisons = list(c("TRhi", "TRlo"))) +
    facet_wrap(~gene_signature, scales="free", nrow=1) +
    scale_fill_manual(values=tr_pal, guide=F) +
    labs(x="Tumor reactivity (per cell)", y="Score") +
    theme_bw()  + theme(legend.position="top")

ggsave("plots/tr_caushi_oliveira_scores.pdf", height=3.5, width=7, useDingbats=F)

# ==== peripheral blood frequency ====
load(file="r_objects/bulk_data.RData")

# only consider blood frequency at scRNA sampling timepoint
tp2_bulk_tcr_bytimepoint <- subset(bulk_tcr_bytimepoint, timepoint=="T2")

# majority phenotype per CD clone
clone_blood_freq <- data.frame(subset(cfi[[]], !is.na(cdr3s_aa)) %>% 
                                   group_by(pat_cdr3, patient, pat_TRB_aa) %>% 
                                   mutate(clonal_tr_score = mean(tumor_reactivity),
                                          clonal_neotcr_cd4 = mean(neotcr_cd4_40),
                                          clonal_neotcr_cd8 = mean(neotcr_cd8_all)
                                         ) %>%
                                   group_by(pat_cdr3, patient, pat_TRB_aa,
                                            pheno_cluster, cd_cluster,
                                            clonal_tr_score, clonal_neotcr_cd4, clonal_neotcr_cd8
                                           ) %>%
                                   summarize(ncells=n(),
                                            ) %>% 
                                   top_n(1, ncells))

# map to blood TCRb data
clone_blood_freq <- dplyr::left_join(clone_blood_freq,
                              tp2_bulk_tcr_bytimepoint[,c("pat_TRB_aa", "sum_freq")],
                              by="pat_TRB_aa")

colnames(clone_blood_freq)[colnames(clone_blood_freq)=="sum_freq"] <- "T2_blood_freq"

# count unmatched blood TCRbeta frequency as 0
clone_blood_freq[is.na(clone_blood_freq)] <- 0

pheno_order_blood_freq <- clone_blood_freq %>% group_by(pheno_cluster) %>% 
    summarize(median_blood_freq=mean(T2_blood_freq)) %>% arrange(median_blood_freq)

# Fig 4A: peripheral frequency per cluster
ggplot(clone_blood_freq, aes(x=factor(pheno_cluster, levels=pheno_order_blood_freq$pheno_cluster), 
                      y=T2_blood_freq, color=pheno_cluster, fill=pheno_cluster)) + 
    geom_boxplot(alpha=0.8, outlier.shape=NA) +
    facet_grid(~cd_cluster, space="free", scales="free") +
    scale_y_log10() +
    scale_size_binned(name="Clone size", breaks=c(1,10,100, 1000, 1500), limits=c(1, 1500)) +
    scale_color_manual(values=cluster_pal, guide=F) +
    scale_fill_manual(values=cluster_pal, guide=F) +
    labs(x="Cluster", y="Blood frequency", title="TCR clonal frequencies in blood (majority phenotype)") +
    theme_classic()  + rotatex

ggsave("plots/blood_frequency_maj_pheno.pdf", height=4, width=7)

# correlation of TR score vs. blood frequency per cluster
tr_vs_pb <- data.frame(clone_blood_freq %>% group_by(patient, pheno_cluster, cd_cluster) %>% 
                       summarize(avg_blood_freq = mean(T2_blood_freq),
                                 avg_tr_score = mean(clonal_tr_score),
                                 avg_neotcr_cd4_score = mean(clonal_neotcr_cd4),
                                 avg_neotcr_cd8_score = mean(clonal_neotcr_cd8)
                                ))

# Fig 4B: CD4 tumor reactivity score vs. peripheral frequency
ggplot(subset(tr_vs_pb, cd_cluster == "CD4"),
       aes(x=avg_neotcr_cd4_score, y=avg_blood_freq)) +
    geom_point(aes(color=pheno_cluster, shape=patient), size=4, alpha=0.7) +
    scale_color_manual(values=cluster_pal, name = "Cluster") +
    scale_y_log10() +
    coord_cartesian(xlim=c(-0.2, 0.2)) +
    geom_smooth(method="lm", alpha=0.1, color='grey', size=1, lty=1) +
    stat_cor(method = "spearman", label.x = -0.05, label.y = log10(5e-05)) +
    labs(x="Mean NeoTCR-CD4 score", y="Mean peripheral blood frequency",
        title = "Majority phenotype per clone") +
    theme_bw() + guides(color=guide_legend(ncol=2))

ggsave("plots/neotcr_cd4_vs_blood_frequency_per_cluster_majpheno.pdf", height=4, width=7)

# Fig 4B: CD8 tumor reactivity score vs. peripheral frequency
ggplot(subset(tr_vs_pb, cd_cluster == "CD8"),
       aes(x=avg_tr_score, y=avg_blood_freq)) +
    geom_point(aes(color=pheno_cluster, shape=patient), size=4, alpha=0.7) +
    scale_color_manual(values=cluster_pal, name = "Cluster") +
    scale_y_log10() +
    coord_cartesian(xlim=c(-0.2, 0.3)) +
    geom_smooth(method="lm", alpha=0.1, color='grey', size=1, lty=1) +
    stat_cor(method = "spearman", label.x = 0.0) +
    labs(x="Mean CD8 tumor-reactivity score", y="Mean peripheral blood frequency",
        title = "Majority phenotype per clone") +
    theme_bw() + guides(color=guide_legend(ncol=2))

ggsave("plots/tr_cd8_vs_blood_frequency_per_cluster_majpheno.pdf", height=4, width=7)

# ==== save data ====
saveRDS(cfi, file="r_objects/cfi.rds")