# # select sulfate reduction pathway orthologs (done in R)
# while read line; do cp Orthogroup_Sequences/$line.fa ../sulf_hammer/$line.faa; done < sulf_red_pathway.csv
#
# # create alignment
# for i in *faa; do clustalo -i $i -o ${i%%.faa}.sto --outfmt=st --threads=4; done
#
# # create hmm profile
# for i in *sto; do hmmbuild ${i%%.sto}.hmm $i; done
#
# # metagenome
#
# datasets summary genome accession PRJNA707313 PRJNA598416 PRJNA362212 PRJNA635214 --annotated --assembly-source GenBank --as-json-lines > metagenome.jsonl
#
# dataformat tsv genome --inputfile metagenome.jsonl --fields accession,annotinfo-featcount-gene-total,assminfo-biosample-bioproject-accession,assminfo-biosample-accession,assminfo-level,assmstats-contig-n50,assmstats-gc-percent,assmstats-number-of-contigs,assmstats-total-sequence-len,organism-name,organism-tax-id,assminfo-biosample-description-comment
#
# datasets download genome accession PRJNA707313 PRJNA598416 PRJNA362212 PRJNA635214 --annotated --assembly-source GenBank --dehydrated --include protein
#
#
# parallel -j 7 hmmscan --tblout {/.}.res -E 1e-3 sulfate_pathway.hmm {} ::: all/*.faa


# useful packages

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

source("~/00_lib/scripts/parse_hmmer.R")

# load metagenome metadata

metagenome_metadata =
  read_tsv("01_inData/metagenome/metagenome_metadata.tsv",
          name_repair = "universal") %>%
  rename("accession" = "Assembly.Accession",
         "n_genes" = "Annotation.Count.Gene.Total",
         "bioproject" = "Assembly.BioSample.BioProject.Accession",
         "n_bases" = "Assembly.Stats.Total.Sequence.Length") %>%
  replace_na(list(bioproject = "PRJNA362212")) %>%
  mutate(
    groups = case_when(bioproject %in% "PRJNA352737" ~ "Pacific Gyre",
                       bioproject = T ~ "Sediments") %>%
      factor(levels = c("Sediments", "Pacific Gyre")),
    bioproject = factor(bioproject, levels = c(
      "PRJNA352737",
      "PRJNA362212",
      "PRJNA598413",
      "PRJNA598416",
      "PRJNA635214",
      "PRJNA707313")),
    n_genome = n(), .by = bioproject)

metagenome_taxonomy =
  taxonomizr::getTaxonomy(pull(metagenome_metadata, Organism.Taxonomic.ID) %>%
                            unique(),
                          "nameNode.sqlite") %>%
  as.data.frame() %>%
  rownames_to_column("Organism.Taxonomic.ID") %>%
  mutate(Organism.Taxonomic.ID = as.double(Organism.Taxonomic.ID))

metagenome_metadata_t =
  left_join(metagenome_metadata, metagenome_taxonomy)

# load experiment results

hmmscan_results =
  fs::dir_ls("01_inData/metagenome/all/hmmsearch_results/",
             glob = "*.tblout") %>%
  map(read_tblout) %>%
  list_rbind(names_to = "accession")


seq_length =
  fs::dir_map("01_inData/metagenome/all/seq_database/",
              fun = \(.x) read_fwf(.x, skip = 1) %>%
                set_names(c("query_name", "seq_length"))) %>%
  list_rbind(names_to = "accession") %>%
  mutate(n_sequences = n()) %>%
  select(-accession)

results_hmmer =
  hmmscan_results %>%
  mutate(accession = str_remove(
    accession,
    "01_inData/metagenome/all/hmmsearch_results/") %>%
           str_remove(".tblout")) %>%
  left_join(metagenome_metadata_t) %>%
  mutate(Ds_variabilis_3be13 = str_remove(domain_name, "-i1"),
         mmx_score = (sequence_score - min(sequence_score)) /
           (max(sequence_score) - min(sequence_score)),
         bioproject2 = recode(bioproject,
                              PRJNA352737 = "Pacific Gyre",
                              PRJNA362212 = "Guaymas Basin",
                              PRJNA598413 = "Gulf of Khambhat",
                              PRJNA598416 = "Gulf of Kutch",
                              PRJNA635214 = "Challenger Deep",
                              PRJNA707313 = "South China Sea") %>%
           as.factor(),
         .by = domain_name) %>%
  left_join(sulf_red_prot, by = c("Ds_variabilis_3be13" = "Locustag")) %>%
  replace_na(list(Gene = "ppaC", Operon = "ppa")) %>%
  left_join(seq_length)

write_csv(results_hmmer, "results_hmmer_uniref.csv")
results_hmmer = read_csv("results_hmmer_uniref.csv")

# for the others

results_tidy =
  results_hmmer %>%
  filter(best_domain_evalue < 1/n_sequences,
         seq_length < (median(seq_length) * 4/3),
         seq_length > (median(seq_length) * 2/3),
         .by = domain_name) %>%
  mutate(presence = 1) %>%
  pivot_wider(id_cols = accession,
              names_from = domain_name,
              values_from = presence,
              values_fn = first,
              values_fill = 0) %>%
  pivot_longer(-accession,
               names_to = "Ds_variabilis_3be13",
               values_to = "presence") %>%
  mutate(Ds_variabilis_3be13 = str_remove(Ds_variabilis_3be13, "-i1")) %>%
  left_join(sulf_red_prot, by = c("Ds_variabilis_3be13" = "Locustag")) %>%
  replace_na(list(Gene = "ppaC", Operon = "ppa")) %>%
  left_join(metagenome_metadata_t)

with(results_tidy, table(order == "Desulfobacterales", presence == 1)) %>%
  prop.test()

res_fisher =
  results_tidy %>%
  nest(data = -Gene) %>%
  mutate(ft = map(data, ~ with(.x,
                               table(order == "Desulfobacterales",
                                     presence == 1)) %>%
                    fisher.test() %>% broom::tidy())) %>%
  unnest(ft)

write_csv(res_fisher %>% select(-data), "res_fisher.csv")

res_filtered =
  results_hmmer %>%
  filter(best_domain_evalue < 1/n_sequences,
         .by = domain_name) %>%
mutate(Operon =
         case_when(
           Op == 8 ~ "*dsr",
           Op == 16 ~ "rnf2",
           Op == 17 ~ "rnf1",
           .default = gsub("_?[A-Z]?[0-9]?", "", Gene)) %>%
         factor(levels = c("nss", "hpp", "ppa", "sat", "apr", "qmo", "dsr", "*dsr", "cdh",
                           "met", "coo", "fol", "fhs", "fdh", "qrc", "rnf2", "rnf1", "nfn",
                           "tmc", "nqr", "bam", "dch", "had", "oah", "bzd")),
       Gene = as.factor(Gene))

res_filtered %>%
  count(domain_name, bioproject2) %>% View()

p1 =
  res_filtered %>%
  group_by(Gene, bioproject2, Operon, groups, n_genome, accession) %>%
  summarise(p = n() > 0) %>%
  summarise(
    p = sum(p),
    n_genome = mean(n_genome),
    prevalence = p/n_genome) %>%
  ggplot() +
  geom_tile(aes(x = bioproject2,
                y = Gene,
                fill = prevalence)) +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        legend.position = "bottom",
        panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0.1, "lines"),
        axis.title.y = element_blank(),
        strip.text.y.right = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(size = 12),
        strip.text = element_text(size = 12)) +
  facet_grid(rows = vars(Operon),
             cols = vars(groups),
             space = "free", scales = "free") +
  scale_fill_gradient(low = "#F5F5F5", high = "#DFC27D")

p1

p2 =
  res_filtered %>%
  group_by(Gene, bioproject2, Operon, groups, order) %>%
  summarise(n = n()) %>%
  summarise(entropy = sum(n/sum(n) * -log(n/sum(n)))) %>%
  ggplot() +
  geom_tile(aes(x = bioproject2,
                y = Gene,
                fill = entropy)) +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        legend.position = "bottom",
        panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0.1, "lines"),
        axis.title.y = element_blank(),
        strip.text.y.right = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 12)) +
  facet_grid(Operon ~ groups,
             space = "free", scales = "free") +
  scale_fill_gradient(low = "#F5F5F5", high = "#80CDC1")


p1 + p2

p3 =
  res_filtered %>%
  mutate(order2 =
           as.factor(order) %>%
           fct_na_value_to_level("Other") %>%
           fct_lump_n(5, other_level = "Other") %>%
           factor(levels = c("Other", "Anaerolineales",
                             "Bacteroidales", "Pseudomonadales",
                             "Desulfobacterales"))) %>%
  ggplot() +
  geom_bar(aes(y = Gene,
               fill = order2),
           position = "stack",
           color = "white") +
theme(axis.text.x = element_text(size = 10, angle = 90),
      legend.position = "bottom",
      panel.spacing.y = unit(0.5, "lines"),
      panel.spacing.x = unit(0.5, "lines"),
      axis.title.y = element_blank(),
      strip.text.y.right = element_blank(),
      axis.title.x = element_blank(),
      panel.border = element_blank(),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 12)) +
  facet_grid(Operon ~ groups + bioproject2, scales = "free",
             space = "free_y") +
  scale_fill_manual(values = c("grey95", "grey80", "grey65", "grey50", "#FFC107")) +
  scale_x_log10()


res_filtered %>%
  mutate(order2 =
           as.factor(order) %>%
           fct_na_value_to_level("Other") %>%
           fct_lump_n(5, other_level = "Other") %>%
           factor(levels = c("Other", "Anaerolineales",
                             "Bacteroidales", "Pseudomonadales",
                             "Desulfobacterales"))) %>%
  filter(bioproject2 %in% c("Guaymas Basin", "Gulf of Kutch",
                            "Challenger Deep", "South China Sea")) %>%
  group_by(order2, bioproject2, Gene, Operon) %>%
  summarise(No_of_detections = n()) %>%
  group_by(order2, Gene, Operon) %>%
  summarise(median = median(No_of_detections),
            sd = sd(No_of_detections)) %>%
  pivot_wider(names_from = order2, values_from = c(median, sd),
              values_fill = 0) %>%
  write_csv("taxonomic_share_median_by_gene.csv")

p3

p1 + p2 + p3 +
  plot_layout(widths = c(1, 1, 1), guides = "collect",axes = "collect_y")

set.seed(1)

future::plan(future::multisession, workers = 6)

my_mod =
  results_tidy %>%
  group_nest(Gene) %>%
  mutate(
    # dwn_data =
    #   furrr::future_map(data,
    #                      ~ group_by(.x, bioproject) %>%
    #                        slice_sample(n = 147),
    #                     .options = furrr::furrr_options(seed = 123)),
    mod =
      furrr::future_map(data,
                         ~ glm(presence ~ bioproject,
                               data = .x,
                               family = "binomial"),
                        .options = furrr::furrr_options(seed = 123)),
    tidies = furrr::future_map(mod, ~ broom::tidy(.x),
                               .options = furrr::furrr_options(seed = 123)))

 my_tidy =
  my_mod %>%
  dplyr::select(Gene, tidies) %>%
  unnest(tidies)

my_data =
  my_mod %>%
  unnest(data) %>%
  group_by(bioproject, Operon, Gene) %>%
  reframe(prevalence = sum(presence > 0)/n_genome,
            presence = sum(presence)) %>%
  distinct()

tidy_2 =
  my_tidy %>%
  pivot_wider(names_from = term, values_from = c(3:6)) %>%
  pivot_longer(cols = -c(1, contains("Intercept")),
               names_to = c(".value", "bioproject"),
               names_sep = "_") %>%
  mutate(bioproject = str_remove(bioproject, "bioproject"),
         padj = p.adjust(p.value, "fdr"),
         bioproject2 = recode(bioproject,
                              PRJNA352737 = "Pacific Gyre",
                              PRJNA362212 = "Guaymas Basin",
                              PRJNA598413 = "Gulf of Khambhat",
                              PRJNA598416 = "Gulf of Kutch",
                              PRJNA635214 = "Challenger Deep",
                              PRJNA707313 = "South China Sea") %>% as.factor(),
         groups = case_when(bioproject %in% "PRJNA352737" ~ "Pacific Gyre",
                            bioproject = T ~ "Sediments") %>%
           factor(levels = c("Sediments", "Pacific Gyre"))) %>%
  left_join(sulf_red_prot) %>%
  left_join(my_data)

pal_met =   c("#FFECB3", "#FFD54F", "#FFC107", "#FFA000", "#FF6F00")


p4 =
  my_tidy %>%
  mutate(intercept = if_else(term %in% "(Intercept)", estimate, NA),
         estimate2 = if_else(!term %in% "(Intercept)", mean(intercept, na.rm = T) + estimate, estimate),
         .by = Gene) %>%
  mutate(bioproject = str_remove(term, "bioproject"),
         padj = p.adjust(p.value, "fdr"),
         bioproject2 = recode(bioproject,
                              `(Intercept)` = "Pacific Gyre",
                              PRJNA362212 = "Guaymas Basin",
                              PRJNA598413 = "Gulf of Khambhat",
                              PRJNA598416 = "Gulf of Kutch",
                              PRJNA635214 = "Challenger Deep",
                              PRJNA707313 = "South China Sea") %>% as.factor(),
         groups = case_when(bioproject2 %in% "Pacific Gyre" ~ "Pacific Gyre",
                            .default = "Sediments") %>%
           factor(levels = c("Pacific Gyre", "Sediments"))) %>%
  left_join(sulf_red_prot) %>%
  ggplot() +
  geom_tile(aes(x = bioproject2,
                y = Gene,
                fill = estimate2)) +
  theme(axis.text.x = element_text(size = 10, angle = 90),
        legend.position = "bottom",
        panel.spacing.y = unit(0, "lines"),
        panel.spacing.x = unit(0.1, "lines"),
        axis.title.y = element_blank(),
        strip.text.y.right = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_blank(),
        strip.text = element_text(size = 12)) +
  facet_grid(rows = vars(Operon),
             cols = vars(groups),
             space = "free", scales = "free") +
  scale_fill_gradient2(high = "#F1A340", mid = "#F7F7F7", low = "#998EC3")

cairo_pdf("example_many_heatmaps.pdf", 24, 16)
p1 + p2 + p4 + p3 +
  plot_layout(guides = "collect", nrow = 1, widths = c(1, 1, 1, 3))
dev.off()

p2 =
  tidy_2 %>%
  filter(!Operon %in% "rfn1") %>%
  ggplot(aes(y = statistic, color = bioproject2)) +
  geom_density(aes(fill = bioproject2, fill = after_scale(colorspace::lighten(fill))),
               alpha = .3) +
  scale_color_manual(values = pal_met,
                     aesthetics = c("color", "fill")) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        panel.grid.major.y = element_line(linewidth = 1, linetype = 3)) +
  scale_y_continuous(breaks = c(10, 5, 0, -5, -10))


cairo_pdf("04_figures/ultimate/metagenome_orange.pdf", 24, 16)
p1 + p2 + plot_layout(widths = c(3, 1))
dev.off()

svg("04_figures/ultimate/metagenome_orange.svg", 24, 16)
p1 + p2 + plot_layout(widths = c(3, 1))
dev.off()



