### Comparative Genomics of 43 Desulfobacteraceae
## Roberto Siani
# 2024

## summary page to use as metadata
#
# datasets summary genome taxon Desulfobacteraceae --assembly-source RefSeq --as-json-lines | dataformat tsv genome --fields accession,annotinfo-featcount-gene-total,assmstats-contig-n50,assmstats-gc-percent,assmstats-number-of-contigs,assmstats-total-sequence-len,organism-name > metadata.tsv
#
## download and rehydrate dataset
#
# datasets download genome taxon Desulfobacteraceae --assembly-source RefSeq --dehydrated
#
# datasets rehydrate --directory ./
#
## move to forever home
#
# cd ncbi_dataset/data && mkdir mkdir ../../fna
#
# for i in GCF*; do mv $i/*fna ../../fna/$i.fna; done
#
# cd ../../fna
#
## remove weird candidatus and dereplicate
#
# Assembly-Dereplicator/dereplicator.py genomes derep
#
# rm GCF_000516475.1.fna
#
# orthofinder -f . -t 5 -og -y -S diamond_ultra_sens -I 3 -a 5

# base packages

pacman::p_load(
  tidyverse,
  patchwork,
  pacman)

# a minimal theme

my_theme <- theme_minimal() +
  theme(plot.margin = margin(3, 5, 0, 3),
        text = element_text(size = 30,
                            family = "Fira Sans",
                            color = "#212121"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(hjust = 1),
        axis.title.y = element_text(hjust = 1),
        legend.position = "bottom",
        legend.justification = c(1, 0),
        legend.margin = margin(0, 5, 3, 5),
        legend.title = element_blank(),
        panel.spacing.x = unit(3, "lines"),
        panel.spacing.y = unit(3, "lines"),
        axis.line = element_blank(),
        panel.border = element_rect(color = "#212121", linewidth = 2 , fill = NA)
  )

theme_set(my_theme)

## load gene counts matrix (MCL I = 3)

geneCounts <-
  read_tsv("data/Orthogroups.GeneCount.tsv") %>%
  mutate(
    prev = rowSums(x = .[, 2:44] > 0)/43,
    partition = case_when(
      prev >= 0.975 ~ "Core",
      prev < 0.975 & prev >= 0.15 ~ "Shell",
      prev < 0.15 & prev > 1/43 ~ "Cloud",
      prev == 1/43 ~ "Unique") %>%
      factor(levels = c("Unique", "Cloud", "Shell", "Core")))

table(geneCounts$partition)

janitor::tabyl(geneCounts, partition)

# load metadata from ncbi datasets

metaData =
  read_tsv("data/metadata.tsv", name_repair = "universal") %>%
  mutate(Organism.Name = gsub("Candidatus|uncultured", "", Organism.Name) %>%
           str_trim(),
         Genome = Assembly.Accession) %>%
  select(Organism.Name, Genome) %>%
  separate(Organism.Name, into = c("Genus", "Species"), sep = " ", remove = F) %>%
  filter(!is.na(Genome)) %>%
  distinct()

metaData[52,] = list(
  "Desulfosarcina variabilis 3be13",
  "Desulfosarcina",
  "variabilis",
  "Ds_variabilis_3be13")

# fix order

metaData = metaData[match(colnames(geneCounts[, 2:44]), metaData$Genome),]

# extract count matrix

panMat = geneCounts[, 1:44] %>% column_to_rownames("Orthogroup") %>% t()

# rarefaction curve and fitting heaps law

rarefaction <- function(pan.matrix, n.perm = 1){
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  nmat <- matrix(0, nrow = nrow(pan.matrix), ncol = n.perm)
  cm <- apply(pan.matrix, 2, cumsum)
  nmat[,1] <- rowSums(cm > 0)
  if(n.perm > 1){
    for(i in 2:n.perm){
      cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)
      nmat[,i] <- rowSums(cm > 0)
      cat(i, "/", n.perm, "\r")
    }
  }
  nmat <- rbind(rep(0, ncol(nmat)), nmat)
  tibble(Genomes = 0:nrow(pan.matrix)) %>%
    bind_cols(as_tibble(nmat, .name_repair = "minimal")) -> rtbl
  colnames(rtbl) <- c("Genome", str_c("Perm", 1:n.perm))
  return(rtbl)
}

heaps <- function(pan.matrix, n.perm = 100){
  objectFun <- function(p, x, y){
    y.hat <- p[1] * x^(-p[2])
    J <- sqrt(sum((y - y.hat)^2))/length(x)
    return(J)
  }
  pan.matrix[which(pan.matrix > 0, arr.ind = T)] <- 1
  ng <- dim(pan.matrix)[1]
  nmat <- matrix(0, nrow = nrow(pan.matrix) - 1, ncol = n.perm)
  for(i in 1:n.perm){
    cm <- apply(pan.matrix[sample(nrow(pan.matrix)),], 2, cumsum)
    nmat[,i] <- rowSums((cm == 1)[2:ng,] & (cm == 0)[1:(ng-1),])
    cat(i, "/", n.perm, "\r")
  }
  x <- rep((2:nrow(pan.matrix)), times = n.perm)
  y <- as.numeric(nmat)
  p0 <- c(mean(y[which(x == 2)] ), 1)
  fit <- optim(p0, objectFun, gr = NULL, x, y, method = "L-BFGS-B", lower = c(0, 0), upper = c(10000, 2))
  p.hat <- fit$par
  names(p.hat) <- c("Intercept", "alpha")
  return(p.hat)
}

# check saturation of the pangenome

set.seed(1)

res_heaps =
  heaps(panMat, n.perm = 999) # open at alpha of ~0.9

set.seed(1)

res_rarefy =
  rarefaction(panMat, 999)

res_rarefy %>%
  pivot_longer(-Genome, names_to = "Permutation", values_to = "Clusters") %>%
  ggplot(aes(x = Genome, y = Clusters)) +
  geom_point(color = "#dddddd") +
  geom_smooth(color = "#BF360C", linewidth = 3) +
  annotate(geom = "text",
           x = 35, y = 100,
           label = paste("alpha =", res_heaps[2] %>% round(2)),
           size = 10,
           color = "#333333") +
  scale_y_continuous(n.breaks = 10) +
  scale_x_continuous(n.breaks = 10) +
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        panel.grid.major.x = element_line(),
        panel.grid.minor.x = element_line())

# TSNE

.Random.seed = readRDS("data/seed_tsne.RDS")

res.tsne =
  Rtsne::Rtsne(dist(panMat, "binary"),
               perplexity = 10,
               initial_dims = 22,
               theta = 0)

factoextra::fviz_nbclust(res.tsne$Y, kmeans)

# clustering

data_sne =
  kmeans(res.tsne$Y, 7) %>%
  broom::augment(res.tsne$Y) %>%
  mutate(Genome = rownames(panMat)) %>%
  left_join(metaData) %>%
  mutate(Genus2 = if_else(str_detect(Genus, "Desulfonema"),
                          str_c(Genus, Species, sep = " "),
                          Genus))

# plot (colors where changed in post-production)

pal_oi = c( "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
            "#E69F00", "#56B4E9","#009E73", "#999999")

plot_3 =
  ggplot(
    data =
      data_sne,
    aes(X1, X2, color = .cluster, label = Genus2, fill = after_scale(color))) +
  scale_shape_manual(values = c(4, 20)) +
  ggrepel::geom_text_repel(size = 12/.pt,
                           aes(fontface = ifelse(str_detect(Genome, "^D"),
                                              "bold",
                                              "plain"))) +
  scale_alpha_manual(values = c(1, .67)) +
  ggalt::geom_encircle(alpha = .2) +
  geom_point(shape = 1, size = 2, stroke = 1) +
  scale_color_manual(values = pal_oi) +
  ggtitle("t-SNE") +
  theme(axis.text = element_blank(), axis.title = element_blank(),
        legend.position = "none") +
  labs(x = "Dim.1", y = "Dim.2")

# subset core genome and match with Ds locustags for interpretability

core_genome =
  read_fwf("data/ds_names.txt",
           fwf_empty("data/ds_names.txt",
                     col_names = c("Ds_variabilis_3be13", "Description"))) %>%
  mutate(Description = str_remove(Description, "Dvar_[0-9]* ")) %>%
  left_join(read_tsv("data/Orthogroups.tsv",
                     col_select = c(1, 6)) %>%
              separate_longer_delim(Ds_variabilis_3be13, ", "),
            multiple = "all") %>%
  select(-Ds_variabilis_3be13) %>%
  left_join(read_tsv("data/Orthogroups.tsv"),
            by = "Orthogroup") %>%
  distinct()

write_tsv(core_genome, "data/core_genes.tsv")

## DCGs genes

dcgs =
  read_tsv("data/dcgs.csv") %>%
  mutate(Operon =
           case_when(
             Op == 8 ~ "*dsr",
             Op == 16 ~ "rnf2",
             .default = gsub("_?[A-Z]?[0-9]?", "", Gene)) %>%
           factor(levels = c("nss", "hpp", "ppa", "sat", "apr", "qmo", "dsr", "*dsr","qrc",
                              "rnf2", "nfn", "tmc", "nqr", "cdh",
                             "met", "coo", "fol", "fhs", "fdh", "pcc", "bm", "mce",
                              "suc", "sdh", "fum", "mae",
                             "mdh", "pyc", "por", "bam", "dch", "had", "oah", "bzd")))

# orthogroups containing DCGs

orthogroups =
  read_tsv("data/Orthogroups.tsv") %>%
  separate_rows(Ds_variabilis_3be13, sep = ", ") %>%
  filter(Ds_variabilis_3be13 %in% dcgs$Locustag |
           grepl("HRM2_16920", `Dt_autotrophicum_HRM2`) |
           grepl("Dmul_07610", `Dc_multivorans_1be1`)) %>%
  left_join(dcgs, by = c("Ds_variabilis_3be13" = "Locustag")) %>%
  mutate(Gene = case_when(`Dt_autotrophicum_HRM2` == "HRM2_16920" ~ "ppaC",
                          `Dc_multivorans_1be1` == "Dmul_07610" ~ "maeA",
                          .default = Gene),
         Description = case_when(`Dt_autotrophicum_HRM2` == "HRM2_16920" ~ "PpaC: Pyrophosphate phospho-hydrolase",
                          `Dc_multivorans_1be1` == "Dmul_07610" ~ "MaeA: NAD-dependent malic enzyme", .default = Description),
         Locustag = case_when(`Dt_autotrophicum_HRM2` == "HRM2_16920" ~ "HRM2_16920",
                              `Dc_multivorans_1be1` == "Dmul_07610" ~ "Dmul_07610",
                              .default = Ds_variabilis_3be13))


write_tsv(orthogroups %>%
            select(Orthogroup, Locustag),
          "data/dcgs_orthogroups",
          col_names = F)

# now select all those orthogroups from Orthogroup_sequences and move them to a separate folder
# to choose only the query DCGs for Uniref, concatenate orthogroups and then write the sequences to separate
# fasta files


target = dcgs$Locustag %>%
  as_vector() %>%
  set_names()


multifasta = seqinr::read.fasta("data/all_dcgs.faa",
                        seqtype = "AA",
                        as.string = T,
                        whole.header = F,
                        strip.desc = F)

# extract sequences

results1 =
  multifasta[target]

map(names(results1),
    ~ seqinr::write.fasta(sequences = results1[.x],
                  names = .x,
                  nbchar = 80,
                  as.string = T,
                  file.out = str_c("data/search_sequences/", .x, ".faa")))

# we build alignment with uniref using jackhmmer
# parallel -j 3 jackhmmer -N 1 -A {.}.sto --incE 1e-9 -E 1e-9 {} ~/00_lib/uniref90.fasta ::: *.faa
# hmmbuild
# hmmpress
# prepare list of orthogroups and list of hmm
# parallel -a ../hmm_id.txt -a ../orthogroups_id.txt --link hmmscan --tblout {2.}.tblout {1} ../dcgs_orthologs/{2}
# from tblout to csv
# for i in *.tblout; do grep -v "^#" $i | awk '{print $1","$3","$5","$6","7","$8","$9","$10}' > ${i%%.tblout}.csv; done

# uniref ------------------------------------------------------------------

# load hmmscan results

ur_tblout =
  list.files("data/search_sequences/",
             pattern = ".csv", full.names = T) %>%
  map(~ read_csv(.x,
                 id = "Orthogroup",
                 col_names = c("target", "query", "evalue", "score"))) %>%
  list_rbind() %>%
  mutate(Orthogroup =
           str_remove(
             Orthogroup,
             "data/search_sequences/") %>%
           str_remove(".csv"),
         Genome = str_remove(query, "_C?[0-9]*$") %>%
           case_match(.,
             "Dvar" ~ "Ds_variabilis_3be13",
             "Dmul" ~ "Dc_multivorans_1be1",
             "dnl" ~ "Dn_limicola",
             "dnm" ~ "Dn_magnum",
             "TOL2" ~ "Db_toluolica_Tol2",
             "HRM2" ~ "Dt_autotrophicum_HRM2",
             .default = .))

write_csv(ur_tblout, "ur_tblout.csv")

ur_tblout = read_csv("ur_tblout.csv")

# prepare for plotting

ur_simple =
  ur_tblout %>%
  group_by(Orthogroup, Genome) %>%
  summarise(score = max(score)) %>%
  reframe(
    Genome,
    score,
    scaled_score = (score - min(score)) / (max(score) - min(score))) %>%
  left_join(orthogroups %>%
              select(Orthogroup, Gene)) %>%
  left_join(metaData %>% select(Genome, Genus, Species)) %>%
  left_join(data_sne) %>%
  left_join(dcgs) %>%
  mutate(.cluster = .cluster %>% fct_rev())

# heatmap plot

ur_simple %>%
  ggplot() +
  geom_point(aes(x = Organism.Name,
                 y = Gene,
                 fill = scaled_score),
             stroke = 1,
             color = "#555555",
             shape = 22,
             size = 7) +
  facet_grid(Operon ~ .cluster,
             space = "free",
             scales = "free", switch = "y") +
  scale_x_discrete(position = "top") +
  scale_fill_binned(type = "gradient",
                    breaks = c(0, 0.25, 0.5, 0.75, 1),
                    low = "seashell",
                    high = "steelblue") +
  theme(panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0, "lines"),
        panel.border = element_rect(linewidth = 0.1),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        strip.text.y.left = element_blank(),
        legend.text = element_text(size = 15),
        legend.key.width = unit(3, "lines"),
        axis.title.y = element_text(hjust = 0, vjust = 1.25, angle = 0),
        axis.text.x.top = element_text(angle = 90, size = 15, hjust = 0),
        axis.text.y = element_text(size = 15, angle = 0),
        strip.placement = "inside")

ggsave("heatmap.svg", width = 12, height = 24)
