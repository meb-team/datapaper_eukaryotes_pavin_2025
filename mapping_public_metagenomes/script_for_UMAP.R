library(umap)
library(compositions)
library(zCompositions)
library(vegan)
library(factoextra)
library(tidyverse)
library(UpSetR)
library(pheatmap)
library(viridisLite)

# 0. Setup ----
my.seed <- 6843651
set.seed(my.seed)
indata <- "data" # to update
outdata <- "results" # to update
if (!dir.exists(outdata)) {
  dir.create(paste0(outdata), recursive = TRUE)
}

# 1. Load metadata ----
### This is an example, use yours here
metadata <- read.csv(paste0(indata, "/metadata_with_reads.tsv"),
  sep = "\t", header = TRUE, row.names = 1, na.strings = ""
)
metadata$country <- as.factor(metadata$country)
metadata$broad_geo_region <- as.factor(metadata$broad_geo_region)
metadata$SAMPLE.MATERIAL <- factor(
  metadata$SAMPLE.MATERIAL,
  levels = c(
    NA, "biofilm (ENVO:000002034)", "brackish water (ENVO:00002019)",
    "brine (ENVO:00003044)", "cryoconite", "estuarine water (ENVO:01000301)",
    "microbial mat (ENVO:01000008)", "sediment (ENVO:00002007)",
    "sediment (ENVO:00002007) biofilm (ENVO:000002034)",
    "sludge (ENVO:00002044)", "stromatolite (ENVO:00002157)",
    "water (ENVO:00002006)", "water (ENVO:00002006) cryoconite",
    "water (ENVO:00002006) sediment (ENVO:00002007)"
  ),
  exclude = NULL
)
metadata$ECOSYSTEM.TYPE <- as.factor(metadata$ECOSYSTEM.TYPE)
metadata$MINIMUM.SIZE.FRACTION <- factor(
  metadata$MINIMUM.SIZE.FRACTION,
  levels = c(
    NA, "0.1 µm", "0.2 µm", "0.22 µm", "0.25 µm", "0.45 µm", "0.65 µm",
    "0.8 µm", "1.2 µm", "2 µm", "2.5 µm", "3 µm", "5 µm", "10 µm",
    "20 µm", "64 µm", "100 µm", "120 µm"
  ), exclude = NULL
)
metadata$MAXIMUM.SIZE.FRACTION <- factor(
  metadata$MAXIMUM.SIZE.FRACTION,
  levels = c(
    NA, "0.2 µm", "0.45 µm", "0.65 µm", "0.7 µm", "0.8 µm", "1 µm",
    "1.2 µm", "1.6 µm", "1.8 µm", "2 µm", "2.5 µm", "3 µm", "4 µm",
    "5 µm", "10 µm", "20 µm", "20.5 µm", "27 µm", "50 µm", "64 µm",
    "100 µm", "120 µm", "156 µm", "180 µm", "297 µm", "~500 µm", "1 mm"
  ),
  exclude = NULL
)
metadata$salinity <- factor(metadata$salinity,
  levels = c(
    NA, "freshwater", "freshwater_sediment",
    "hot_spring", "alkaline_water",
    "estuarine_water", "estuarine_sediment",
    "saline_water", "saline_sediment"
  )
)

metadata$clean_reads <- as.numeric(metadata$clean_reads)
metadata$clean_bases <- as.numeric(metadata$clean_bases)
metadata$mapped_reads <- as.numeric(metadata$mapped_reads)
metadata$filtered_reads <- as.numeric(metadata$filtered_reads)

# 2. Load number of reads ----
df_reads <- read.table(
  paste(
    indata,
    "sum_read_per_MAG.tsv",
    sep = "/"
  ),
  sep = "\t", row.names = 1, header = TRUE
)

# Check the dataframe
df_reads[1:10, c(1:5, 452:455)]

## Remove MAGs
### A command to identify the number of read per MAG :
###     df_reads %>% rowSums() %>% sort() %>% head()
### Two have 1500 or less reads mapped, they add nothing in the data, drop them
mag_to_drop <- to_delete <- rownames(df_reads)[rowSums(df_reads) < 1500]
df_reads <- df_reads[row.names(df_reads) != mag_to_drop, ]

## Remove metagenomes
### When =< 100 reads mapped
metag_to_drop <- colnames(df_reads)[colSums(df_reads) < 100]
df_reads <- df_reads[, !colnames(df_reads) %in% metag_to_drop]

# 3. Detection ----
### Import the data about MAG detection and keep the same rows and columns as
### in the table above (number of reads)
df_detec <- read.table(paste0(indata, "/detection_per_MAG.tsv"),
  sep = "\t", row.names = 1, header = TRUE
)

df_detec <- df_detec[rownames(df_reads), colnames(df_reads)]

# 4. Import taxonomy ----
## Headers = id	superkingdom	clade	phylum	class	order	family	genus	species
df_taxo <- read.table(paste0(indata, "/taxo_MAG.tsv"),
  sep = "\t", row.names = NULL, header = TRUE, na.strings = ""
)

df_taxo_phylum <- df_taxo[, c("id", "phylum")]

# 5. The UMAP ----
## Do the UMAP -----
df_reads_clr <- as.data.frame(
  clr(cmultRepl(df_reads, method = "CZM", z.warning = 1))
)

custom_config <- umap.defaults
custom_config$n_neighbors <- 25
custom_config$min_dist <- 0.25
# default 15 and 0.1 respectivelly

df_reads_clr_umap <- umap(t(df_reads_clr), config = custom_config)

## Save/load the UMAP -----
saveRDS(df_reads_clr_umap, paste0(outdata, "/", "UMAP_nb_reads_CLR.RDS"))
# df_reads_clr_umap <- readRDS(paste0(outdata, "/", "UMAP_nb_reads_CLR.RDS"))

## Get the coordinates in the UMAP space
df_reads_clr_umap_df <- as.data.frame(df_reads_clr_umap$layout)
colnames(df_reads_clr_umap_df) <- c("umap1", "umap2")

## Add metadata in the dataframe
df_reads_clr_umap_df <- merge(df_reads_clr_umap_df, metadata,
  by = "row.names",
  all.x = TRUE
)
rownames(df_reads_clr_umap_df) <- df_reads_clr_umap_df$Row.names
df_reads_clr_umap_df$Row.names <- NULL
head(df_reads_clr_umap_df)

# 6. Determine the number of cluster in th UMAP
## Sources :
##  - https://uc-r.github.io/hc_clustering
##  - https://uc-r.github.io/kmeans_clustering
##  - https://www.datanovia.com/en/lessons/determining-the-optimal-number-of-clusters-3-must-know-methods/

## function "hcut" uses "Ward.D2" and Euclidean distance by default
## Elbow methob
fviz_nbclust(df_reads_clr_umap$data, FUN = hcut, method = "wss", k.max = 10)

## Average silhouette method
fviz_nbclust(df_reads_clr_umap$data, FUN = hcut, method = "silhouette", k.max = 10)

## Identify the clusters
## Aitchison distance (euclidean on CLR); HC with Ward.D2 and cut the tree
df_reads_clr_dist <- vegdist(df_reads_clr_umap$data, method = "euclidean")
df_reads_clr_clust <- hclust(df_reads_clr_dist, method = "ward.D2")
clusters_df$three_clusters <- cutree(df_reads_clr_clust, k = 3) # 3 clusters
clusters_df$three_clusters <- as.factor(clusters_df$three_clusters)
clusters_df$metagenome <- rownames(clusters_df)
head(clusters_df)

## Add the information in the dataframe
df_reads_clr_umap_df <- df_reads_clr_umap_df %>%
  rownames_to_column("metagenome") %>%
  as_tibble() %>%
  inner_join(., clusters_df, by = "metagenome") %>% # merge both DF
  mutate(three_clusters = paste0("cluster_", three_clusters))

# 7.1. UpSet plot with shared detected MAG per cluster UMAP cluster ----
# For this I need the number of detected MAG per metagenome

# For the MAG (rows) and clusters of metagenomes (columns), number of metagenome
## from cluster X in which the MAG is found detected at >= 10%

df_detec_per_cluster <- as.data.frame(ifelse(df_detec >= 10, 1, 0))
df_detec_per_cluster <- df_detec_per_cluster %>%
  rownames_to_column("MAG") %>%
  # Make a looooonggg list
  pivot_longer(-MAG, names_to = "metagenome", values_to = "is_MAG_detected") %>%
  left_join(., clusters_df, by = "metagenome") %>% # merge cluster ID
  mutate(three_clusters = paste0("cluster_", three_clusters)) %>%
  group_by(MAG, three_clusters) %>% # Group MAG + cluster
  # Count in how many Metagenome in each cluster a MAG is present
  summarise(present = sum(is_MAG_detected == 1)) %>%
  # Clean table
  pivot_wider(names_from = "three_clusters", values_from = "present") %>%
  as.data.frame()

## Then the UpSet plot
rownames(df_detec_per_cluster) <- df_detec_per_cluster$MAG
df_detec_per_cluster$MAG <- NULL
df_detec_per_cluster_upset <- as.data.frame(
  ifelse(df_detec_per_cluster > 0, 1, 0)
)

### Uncomment the lines bellow to save the plot as SVG
# svg(file = paste0(outdata, "/", "upSet_number_MAG_detected_per_UMAP_cluster.svg"),
#     width = 6, height = 4.5)
upset(df_detec_per_cluster_upset,
  mainbar.y.label = "No. MAGs shared",
  sets.x.label = "No. MAG found"
)
# dev.off()

# 7.2 UpSet Plot with taxonomy ----
### Add taxonomy to the "df_detec_per_cluster" df
df_detec_per_cluster_taxo <- as.data.frame(ifelse(df_detec >= 10, 1, 0))
df_detec_per_cluster_taxo <- df_detec_per_cluster_taxo %>%
  rownames_to_column(var = "id") %>%
  left_join(., df_taxo_phylum, by = "id") %>% # Add taxonomy
  pivot_longer(-c(id, phylum),
    names_to = "metagenome",
    values_to = "is_MAG_detected"
  ) %>% # Pivot longer
  left_join(., clusters_df, by = "metagenome") %>% # add cluster_id
  mutate(three_clusters = paste0("cluster_", three_clusters)) %>%
  group_by(phylum, three_clusters) %>% # Group by phylum + cluster
  # Count in how many Metagenome in each cluster a phylum is present
  summarise(present = sum(is_MAG_detected == 1)) %>%
  # Clean table
  pivot_wider(names_from = "three_clusters", values_from = "present") %>%
  as.data.frame()

## Then the UpSet plot
# replace the "NA"
df_detec_per_cluster_taxo[is.na(df_detec_per_cluster_taxo$phylum), "phylum"] <- "unknown"
# Move a col as rownames and make it binary
rownames(df_detec_per_cluster_taxo) <- df_detec_per_cluster_taxo$phylum
df_detec_per_cluster_taxo$phylum <- NULL
df_detec_per_cluster_taxo_upset <- as.data.frame(
  ifelse(df_detec_per_cluster_taxo > 0, 1, 0)
)

### Uncomment the lines bellow to save the plot as SVG
# svg(file = paste0(outdata, "/", "upSet_phyla_detected_per_UMAP_cluster.svg"),
#    width = 6, height = 4.5)
upset(df_detec_per_cluster_taxo_upset,
  mainbar.y.label = "No. phyla shared",
  sets.x.label = "No. phyla found"
)
# dev.off()

# 7.3 Contingency table ----
## Cluster / vs Environment ; Number of Metagenomes in which at least a MAG is detected

contingency_table <- as.data.frame(ifelse(df_detec >= 10, 1, 0))
contingency_table <- contingency_table %>%
  rownames_to_column(var = "MAG") %>%
  # Make a looooonggg list
  pivot_longer(-MAG, names_to = "metagenome", values_to = "is_MAG_detected") %>%
  left_join(., clusters_df, by = "metagenome") %>% # merge cluster ID
  mutate(three_clusters = paste0("cluster_", three_clusters)) %>%
  left_join(.,
    rownames_to_column(metadata[, c("origin_simple", "ECOSYSTEM.TYPE")]),
    by = join_by("metagenome" == "rowname")
  ) %>% # merge cluster ID
  group_by(three_clusters, ECOSYSTEM.TYPE) %>% # Group MAG + cluster
  distinct(metagenome, is_MAG_detected) %>% # drop duplicates ==> Number of
  # metagenome in which a MAG is detected. NOT the SUM of MAG!
  summarise(present = sum(is_MAG_detected == 1)) %>%
  pivot_wider(
    names_from = "ECOSYSTEM.TYPE", values_from = "present",
    values_fill = 0
  ) %>% # Clean table
  as.data.frame()
rownames(contingency_table) <- contingency_table$three_clusters
contingency_table$three_clusters <- NULL
contingency_table <- as.matrix(contingency_table)

### Uncomment the lines bellow to save the plot as SVG
# svg(file = paste0(outdata, "/",
#   "contingency_nbr_metagenome_with_MAG_detected_par_UMAP_cluster.svg"
# ), width = 10, height = 6)
pheatmap(contingency_table,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  scale = "none",
  number_format = "%i",
  color = viridis(250),
  number_color = "darkgrey",
  cellwidth = 25,
  cellheight = 25,
  fontsize = 12
)
# dev.off()

## Contingency table of the number of metagenome in each category per cluster
metagenome_per_clust_env <- df_reads_clr_umap_df %>%
  select(metagenome, ECOSYSTEM.TYPE, three_clusters) %>%
  group_by(three_clusters, ECOSYSTEM.TYPE) %>%
  summarise(nb_metagenome = n()) %>%
  pivot_wider(
    names_from = ECOSYSTEM.TYPE, values_from = nb_metagenome,
    values_fill = 0
  ) %>%
  as.data.frame()
rownames(metagenome_per_clust_env) <- metagenome_per_clust_env$three_clusters
metagenome_per_clust_env$three_clusters <- NULL

# svg(file = paste0(outdata, "/", "contingency_all_metagenomes.svg"),
#     width = 10, height = 6)
pheatmap(metagenome_per_clust_env,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  scale = "none",
  number_format = "%i",
  color = viridis(250),
  number_color = "darkgrey",
  cellwidth = 25,
  cellheight = 25,
  fontsize = 12
)
# dev.off()

# 10. Export some information ----
## Metagenomes UMAP corrdinates; metadata; cluster ID; number MAG detected
# write.table(df_reads_clr_umap_df,
#             paste0(outfinal, "/metagenomes_UMAP_and_clusters.tsv"),
#             sep = "\t", row.names = FALSE, quote = FALSE)
