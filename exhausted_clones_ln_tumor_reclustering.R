# ==== set up ====
packages <- c("Seurat", "magrittr", "dplyr", "reshape2", "pheatmap", "tidyr",
              "ggplot2", "ggrepel", "ggsci", "ggpubr", "ggrastr",
              "RColorBrewer", "wesanderson", "viridis",
              "gridExtra", "grid", "cowplot")

invisible(lapply(packages, library, character.only = TRUE))

source("~/oak/scripts/sc_utils.R")

# load data
cfi <- readRDS(file="r_objects/cfi.rds")
load(file="r_objects/all_color_pals.RData")

# ==== recluster exhausted tumor-LN CD8 clones ====
# identify exhausted CD8 clones
cd8_tumor_clone_df <- readRDS(file="r_objects/cd8_tumor_clone_df.rds")

dys_hi_clones <- subset(cd8_tumor_clone_df, dys_status == "dys_hi")$pat_cdr3
cfi$dys_hi <- cfi$pat_cdr3 %in% dys_hi_clones

clone_df <- subset(cfi[[]], dys_hi & cd_cluster == "CD8") %>% 
    group_by(patient, pat_cdr3) %>% 
    summarize(clone_size=n(), 
              region_types=list(unique(region_type)),
              nregion_types=length(unique(region_type)),
              ln_present="LN" %in% unique(region_type),
              tumor_present="Tumor" %in% unique(region_type)
             ) %>%
 filter(nregion_types >1)

head(clone_df)

# subset to expanded exhausted CD8 tumor clones found in LN
ln_tumor_expanded_clones <- subset(clone_df, ln_present & tumor_present & clone_size > 2)

# subset scRNA/TCR data to cells in exhausted tumor-LN clones
ln_tumor_ex_cells <- subset(cfi, pat_cdr3 %in% ln_tumor_expanded_clones$pat_cdr3 & region_type != "Adrenal")
ln_tumor_ex_cells

# split cells by patient/batch
ln_tumor_ex_cells$batch <- ln_tumor_ex_cells$patient
ln_tumor_ex_cells$batch[ln_tumor_ex_cells$orig.ident %in% c("MSK1263_R7", "MSK1263_R8")] <- "MSK1263_batch2"

obj_list <- SplitObject(ln_tumor_ex_cells, split.by = "batch")

# integrate data
obj_list <- lapply(X = obj_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
    VariableFeatures(x) <- remove_genes(VariableFeatures(x), species = "human")
    x
})

features <- SelectIntegrationFeatures(object.list = obj_list)
anchors <- FindIntegrationAnchors(object.list = obj_list, anchor.features = features)
ln_tumor_ex_int <- IntegrateData(anchorset = anchors)

# rescale and cluster
DefaultAssay(ln_tumor_ex_int) <- "integrated"

ln_tumor_ex_int <- ScaleData(ln_tumor_ex_int, verbose = FALSE) %>%
    RunPCA(npcs = 40, verbose = FALSE)

ElbowPlot(ln_tumor_ex_int, ndims = 40)

ln_tumor_ex_int <- RunUMAP(ln_tumor_ex_int, reduction = "pca", dims = 1:30) %>%
    FindNeighbors(reduction = "pca", dims = 1:30) %>%
    FindClusters(resolution = 0.3)

pal1 <- pal_locuszoom()(7)[c(1,2,4,7,5,6,3)]

DimPlot(ln_tumor_ex_int, label=T, cols = pal1)
DimPlot(ln_tumor_ex_int, split.by="region_type", cols = pal1)

# find cluster markers
features_to_test <- remove_genes(rownames(ln_tumor_ex_int@assays$RNA), species = "human")

cluster_markers <- FindAllMarkers(ln_tumor_ex_int, assay = "RNA", only.pos = T, 
                                  features = features_to_test)

saveRDS(cluster_markers, file="r_objects/ln_tumor_ex_int_cluster_markers.rds")

top_genes <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt=avg_logFC)

Idents(ln_tumor_ex_int) <- "seurat_clusters"
cluster_averages <- AverageExpression(ln_tumor_ex_int, assays = "RNA", return.seurat = T)

pheatmap(cluster_averages@assays$RNA[unique(top_genes$gene),], angle=0, scale="row", 
         cluster_cols = F, cluster_rows = F,
         treeheight_row = 15, treeheight_col = 15)

DefaultAssay(ln_tumor_ex_int) <- "RNA"
FeaturePlot(ln_tumor_ex_int, 
            features = c("TCF7", "SELL", "IL7R", "CCR7", "ZEB2", "TOX", "PDCD1", "LAG3"),
            ncol = 4)

FeaturePlot(ln_tumor_ex_int,
            features = c("GZMK", "GZMA", "IFITM1", "KLF12", "CD27", "IL32", "NKG7", "RARRES3"), 
            ncol = 4)

# define marker sets
clust4_genes <- subset(cluster_markers, cluster==4 & p_val_adj < 0.05)$gene
exhaustion_genes <- c("TOX", "PDCD1", "HAVCR2", "TIGIT", "ENTPD1", "LAG3", "ID2")
memory_genes <- c("TCF7", "SELL", "IL7R",  "EOMES", "CCR7", "LEF1", "ZNF683")
effector_genes <- c("KLRG1", "KLRC1", "GZMA", "GZMB", "IFNG")
genes_to_plot <- c(clust4_genes, "MKI67", "CX3CR1", "S1PR5", memory_genes, exhaustion_genes)
other_genes <- c("CD28", "CD226", "MKI67", "GZMB", "PRF1", "IFNG")

# cluster annotation
cluster_annotation <- c(`0`="3: TRM", `1`="4: Toxhigh exhausted", `2`="6: TIM-3high exhausted", 
                        `3`="7: Proliferating exhausted", `4`="2: Progenitor exhausted", 
                        `5`="5: PD-1high exhausted", `6`="1: Central memory")

ln_tumor_ex_int$recluster_pheno <- factor(cluster_annotation[ln_tumor_ex_int$seurat_clusters],
                                          levels=unname(sort(cluster_annotation)))

recluster_pheno_pal <- pal_locuszoom()(7)[c(4,2,1,7,5,6,3)]
names(recluster_pheno_pal) <- cluster_annotation

recluster_pheno_pal <- recluster_pheno_pal[order(names(recluster_pheno_pal))]
recluster_pheno_pal['5: PD-1high exhausted'] = "#ff858d"
recluster_pheno_pal['7: Proliferating exhausted']  = "#a83262"

DimPlot(ln_tumor_ex_int, label=F, group.by="recluster_pheno", cols = recluster_pheno_pal)

all_cells <- FetchData(ln_tumor_ex_int,
                      vars=c("UMAP_1", "UMAP_2", "recluster_pheno", "region_type"))

bg_cells <- FetchData(SubsetData(ln_tumor_ex_int, max.cells.per.ident = 3000, random.seed=10), 
                      vars=c("UMAP_1", "UMAP_2", "recluster_pheno", "region_type"))

# Fig 2A: reclustered UMAP
ggplot(all_cells, aes(x=UMAP_1, y=UMAP_2)) +
    rasterise(geom_point(aes(color=recluster_pheno), size=1), dpi=300) +
    scale_color_manual(values=recluster_pheno_pal) +
    theme_classic() + labs(title="All cells")

ggsave("plots/recluster_umap_newcolor.pdf", width=6.2, height=4)

# Fig 2B (top): reclustered UMAP split by region type
cur_rt = "LN"
ggplot(bg_cells, aes(x=UMAP_1, y=UMAP_2)) +
    rasterise(geom_point(color='grey90', size=1), dpi=300) +
    rasterise(geom_point(data=subset(all_cells, region_type == cur_rt), 
               aes(color=recluster_pheno), size=1.5), dpi=300) +
    scale_color_manual(values=recluster_pheno_pal) +
    theme_classic() + labs(title=cur_rt)

ggsave("plots/recluster_ln_umap_newcolor.pdf", width=6.2, height=4)

cur_rt = "Tumor"
ggplot(bg_cells, aes(x=UMAP_1, y=UMAP_2)) +
    rasterise(geom_point(color='grey90', size=1), dpi=300) +
    rasterise(geom_point(data=subset(all_cells, region_type == cur_rt), 
               aes(color=recluster_pheno), size=1.5), dpi=300) +
    scale_color_manual(values=recluster_pheno_pal) +
    theme_classic() + labs(title=cur_rt)

ggsave("plots/recluster_tumor_umap_newcolor.pdf", width=6.2, height=4)

# Fig 2C: cluster gene expression heatmaps
Idents(ln_tumor_ex_int) <- "recluster_pheno"
cluster_averages_annotated <- AverageExpression(ln_tumor_ex_int, assays = "RNA", return.seurat = T)

pheatmap(t(as.matrix(cluster_averages_annotated@assays$RNA[exhaustion_genes,])), cluster_rows = F,
        color = rev(colorRampPalette(brewer.pal(n = 10, "RdYlBu"))(20)),
         angle=45, scale="column",
         treeheight_col = 15,
         filename="plots/recluster_heatmap_exgenes.pdf", width=4.5, height=3
         )

pheatmap(t(as.matrix(cluster_averages_annotated@assays$RNA[memory_genes,])), cluster_rows = F, 
         color = rev(colorRampPalette(brewer.pal(n = 10, "RdYlBu"))(20)),
         angle=45, scale="column",
         treeheight_col = 15,
         filename="plots/recluster_heatmap_memgenes.pdf", width=4.5, height=3
        )

pheatmap(t(as.matrix(cluster_averages_annotated@assays$RNA[other_genes,])), cluster_rows = F, 
         color = rev(colorRampPalette(brewer.pal(n = 10, "RdYlBu"))(20)),
         angle=45, scale="column",
         treeheight_col = 15,
         filename="plots/recluster_heatmap_othergenes.pdf", width=4, height=3
        )

pheatmap(t(as.matrix(cluster_averages_annotated@assays$RNA[clust4_genes,])), cluster_rows = F,
         angle=45, scale="column",
         color = rev(colorRampPalette(brewer.pal(n = 10, "RdYlBu"))(20)),
         treeheight_col = 15,
         filename="plots/recluster_heatmap_clust4genes.pdf", width=13, height=3
        )

pheatmap(cluster_averages_annotated@assays$RNA[effector_genes,], cluster_cols = F,
         angle=45, scale="row",
         treeheight_col = 15
        )

pheatmap(cluster_averages_annotated@assays$RNA[genes_to_plot,], angle=45, scale="row",
         treeheight_row = 15, treeheight_col = 15
        )

table(ln_tumor_ex_int[[c("pheno_cluster", "seurat_clusters")]])
prop.table(table(ln_tumor_ex_int[[c("pheno_cluster", "seurat_clusters")]]), 2) * 100

# cluster proportions per region type
ggplot(ln_tumor_ex_int[[]], aes(x=region_type)) +
    geom_bar(aes(fill=recluster_pheno), position = "fill") +
    scale_fill_manual(values=pal2, name="Cluster") +
    coord_cartesian(expand=F) +
    labs(x="Region type", y="Proportion of cells") +
    theme_bw() + nobg + rotatex

ggsave("plots/recluster_region_pheno_bar.pdf", width=4.5, height=4)

# Fig 2B (bottom): reclustered phenotype proportion bar plot
ggplot(subset(ln_tumor_ex_int[[]], region_type != 'Normal'), aes(x=region_type)) +
    geom_bar(aes(fill=recluster_pheno), position = "fill") +
    scale_fill_manual(values=recluster_pheno_pal, name="Cluster") +
    coord_flip(expand=F) +
    labs(x="Region type", y="Proportion of cells") +
    theme_bw() + nobg #+ rotatex

ggsave("plots/recluster_region_pheno_bar_horizontal_newcolor.pdf", width=8, height=2)

# LN tumor phenotype bar plots
cd8_dys_clones_wln <- readRDS(file="r_objects/cd8_dys_clones_wln.rds")

cd8_dys_clones_wln_noadrenal <- data.frame(subset(cd8_dys_clones_wln, region_type != "Adrenal") %>%
                                           group_by(pat_cdr3) %>% 
                                           mutate(clone_size_noadrenal = sum(ncells)) %>%
                                           arrange(desc(clone_size_noadrenal))
                                           )
head(cd8_dys_clones_wln_noadrenal, 10)

# find top expanded clones
top_clones <- unique(cd8_dys_clones_wln_noadrenal[,c("pat_cdr3","patient","region_types","clone_size_noadrenal")]) %>% 
    group_by(patient) %>% 
    mutate(clone_name=paste0(patient, " clone", rank(desc(clone_size_noadrenal), ties.method = "random"))) %>% 
    top_n(n = 40, wt=clone_size_noadrenal)

top_clones$clone_name <- factor(top_clones$clone_name, 
                                levels=top_clones$clone_name) # rename clones for better plotting

# cells from top expanded clones
dys_clone_cells <- subset(ln_tumor_ex_int, pat_cdr3 %in% unique(top_clones$pat_cdr3))
dys_clone_cells$region_broad <- as.character(dys_clone_cells$region_type)

df <- subset(dys_clone_cells[[]], 
             region_type != "Normal" &
             cd_cluster == "CD8" & pat_cdr3 %in% top_clones$pat_cdr3[1:20]) %>% 
    group_by(pat_cdr3, patient, region_broad) %>% mutate(ncells=n()) 
df <- dplyr::left_join(df, top_clones[,c("pat_cdr3", "clone_name")], by="pat_cdr3")

df$pheno_cluster <- factor(df$pheno_cluster, levels=levels(cfi$pheno_cluster))

# Fig S7D: clonal phenotype bar plots of LN and tumor cells
ggplot(df,
       aes(x=region_broad)) +
    geom_bar(aes(fill=pheno_cluster), position="fill") +
    geom_text(data=unique(df[,c("clone_name", "patient", "region_broad", "ncells")]), 
              aes(label=ncells), y=1.05, size=2) +
    facet_wrap(~clone_name, nrow=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) +
    scale_fill_manual(values=cluster_pal, name="Cluster", limits=force) +
    coord_cartesian(ylim=c(0,1.1), expand=F) +
    labs(x="Region type", y="Proportion of cells in region", 
         title="Top 20 expanded exhausted CD8 clones (found in tumor and LN)") +
    theme_bw() + rotatex + nobg + guides(fill=guide_legend(nrow=1)) + theme(legend.position="bottom")

ggsave(file="plots/ln_tumor_clone_cluster_barplot_broad_region_noadrenal.pdf", width=16, height=3.5)

ggplot(df,
       aes(x=region_broad)) +
    geom_bar(aes(fill=recluster_pheno), position="fill") +
    geom_text(data=unique(df[,c("clone_name", "patient", "region_broad", "ncells")]), 
              aes(label=ncells), y=1.05, size=2) +
    facet_wrap(~clone_name, nrow=1, labeller = label_wrap_gen(width = 2, multi_line = TRUE)) +
    scale_fill_manual(values=recluster_pheno_pal, name="Cluster", limits=force) +
    coord_cartesian(ylim=c(0,1.1), expand=F) +
    labs(x="Region type", y="Proportion of cells in region", 
         title="Top 20 expanded exhausted CD8 clones (found in tumor and LN)") +
    theme_bw() + rotatex + nobg + guides(fill=guide_legend(nrow=1)) + theme(legend.position="bottom")

ggsave(file="plots/ln_tumor_clone_cluster_barplot_recluster_noadrenal.pdf", width=16, height=3.5)

# clonal proportion of LN/tumor cells in progenitor cluster
table(ln_tumor_ex_int[[c("recluster_pheno", "region_type")]])

prop.table(table(ln_tumor_ex_int[[c("recluster_pheno", "region_type")]]), margin=1)*100
prop.table(table(ln_tumor_ex_int[[c("recluster_pheno", "region_type")]]), margin=2)*100

cluster_breakdown <- data.frame(ln_tumor_ex_int[[]] %>% 
                                group_by(patient, pat_cdr3, clone_size, region_type, recluster_pheno, .drop=F) %>%
                                summarize(ncells=n())
                               )

cluster_breakdown$prop <- (cluster_breakdown$ncells / cluster_breakdown$clone_size)*100
head(cluster_breakdown)

ln_progen_clones <- subset(cluster_breakdown, 
                           region_type == "LN" & recluster_pheno == "2: Progenitor exhausted" & ncells != 0)
nrow(ln_progen_clones)
n_distinct(cluster_breakdown$pat_cdr3)
nrow(ln_progen_clones)/n_distinct(cluster_breakdown$pat_cdr3)
table(ln_progen_clones$patient)
table(unique(ln_tumor_ex_int[[c("patient", "pat_cdr3")]])$patient)
table(ln_progen_clones$patient)/table(unique(ln_tumor_ex_int[[c("patient", "pat_cdr3")]])$patient)

df <- data.frame(patient=rep(c("MSK1263", "MSK1302"), each=2), 
                 nclones=c(40, 90-40, 23, 25-23), 
                 category=rep(c("in_progen", "notin_progen"), 2))

sum(df$nclones)
df %>% group_by(patient) %>% summarize(sum(nclones))

# Fig S7B: proportion of clones with LN progenitor cell
cur_pat <- "MSK1263"
p1 <- ggplot(subset(df, patient==cur_pat), 
       aes(x="", y=nclones, fill=category)) +
  geom_bar(stat="identity", width=1) +
  geom_text(aes(label= nclones)) +
  coord_polar("y", start=0, direction = -1) +
  theme_void() + ggtitle(cur_pat)

cur_pat <- "MSK1302"
p2 <- ggplot(subset(df, patient==cur_pat), 
       aes(x="", y=nclones, fill=category)) +
    geom_bar(stat="identity", width=1) +
    geom_text(aes(label=nclones)) +
    coord_polar("y", start=0, direction=-1) +
    theme_void() + ggtitle(cur_pat)

pdf("plots/recluster_ln_progen_pie.pdf", width=6, height=3)
grid.arrange(p1, p2, nrow=1)
dev.off()

# Fig S7C: % of clone in progenitor cluster
ggplot(subset(cluster_breakdown, region_type != "Normal" & recluster_pheno == "2: Progenitor exhausted" ), 
       aes(x=prop, fill=region_type, color=region_type)) +
    stat_bin(alpha=0.7) +
    facet_wrap(~region_type, nrow=2, scales="free_y") +
    scale_color_manual(values=region_type_pal, name="Region type", limits=force) +
    scale_fill_manual(values=region_type_pal, name="Region type", limits=force) +
    labs(x="% of clone in cluster 2: Progenitor exhausted", y="# Clones") +
    coord_cartesian(xlim=c(0,80)) +
    theme_bw()

ggsave("plots/recluster_pct_in_ln_progen.pdf", width=6, height=4)

summary(subset(cluster_breakdown, region_type == "LN" & recluster_pheno == "2: Progenitor exhausted")$prop)
summary(subset(cluster_breakdown, region_type == "Tumor" & recluster_pheno == "2: Progenitor exhausted")$prop)

# save data
saveRDS(ln_tumor_ex_int, file="r_objects/ln_tumor_ex_int.rds")