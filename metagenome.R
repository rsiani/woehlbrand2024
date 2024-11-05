### Quantitative and taxonomic distribution of DCGs in sediment metagenomes
## Roberto Siani
# 2024

# datasets summary genome accession PRJNA707313 PRJNA598416 PRJNA362212 PRJNA635214 --annotated --assembly-source GenBank --as-json-lines > metagenome.jsonl

# dataformat tsv genome --inputfile metagenome.jsonl --fields accession,annotinfo-featcount-gene-total,assminfo-biosample-bioproject-accession,assminfo-biosample-accession,assminfo-level,assmstats-contig-n50,assmstats-gc-percent,assmstats-number-of-contigs,assmstats-total-sequence-len,organism-name,organism-tax-id,assminfo-biosample-description-comment

# datasets download genome accession PRJNA707313 PRJNA598416 PRJNA362212 PRJNA635214 --annotated --assembly-source GenBank --dehydrated --include protein


# parallel -j 3 hmmsearch --tblout {.}.tblout ../dcg.hmm {} ::: *.faa
# for i in *.tblout; do grep -v "^#" $i | awk '{print $1","$3","$5","$6","7","$8","$9","$10}' > ${i%%.tblout}.csv; done

# IMG dataset was manually downloaded from JGI IMG/MER, along with metadata

# read NCBI MAGs metadata

metadata_ncbi =
  read_tsv("data/metadata_ncbi.tsv",
          name_repair = "universal") %>%
  select("accession" = "Assembly.Accession",
         "n_genes" = "Annotation.Count.Gene.Total",
         "bioproject" = "Assembly.BioSample.BioProject.Accession",
         "n_bases" = "Assembly.Stats.Total.Sequence.Length",
         Organism.Taxonomic.ID) %>%
  replace_na(list(bioproject = "PRJNA362212")) %>%
  mutate(
    bioproject = factor(bioproject, levels = c(
      "PRJNA352737",
      "PRJNA362212",
      "PRJNA598413",
      "PRJNA598416",
      "PRJNA635214",
      "PRJNA707313")),
    site = recode(bioproject,
                      PRJNA352737 = "Pacific Gyre",
                      PRJNA362212 = "Guaymas Basin",
                      PRJNA598413 = "Gulf of Khambhat",
                      PRJNA598416 = "Gulf of Kutch",
                      PRJNA635214 = "Challenger Deep",
                      PRJNA707313 = "South China Sea") %>%
      as.factor(),
    n_genome = n(), .by = bioproject)

# taxonomic assignment based on the taxa ID from NCBI

metagenome_taxonomy =
  taxonomizr::getTaxonomy(pull(metadata_ncbi, Organism.Taxonomic.ID) %>%
                            unique(),
                          "nameNode.sqlite") %>%
  as.data.frame() %>%
  rownames_to_column("Organism.Taxonomic.ID") %>%
  mutate(Organism.Taxonomic.ID = as.double(Organism.Taxonomic.ID))

# join

metagenome_ncbi =
  left_join(metadata_ncbi, metagenome_taxonomy) %>%
  select(-Organism.Taxonomic.ID)

# read IMG MAGs metadata
# we add information on the project contained in separate file
# and reformat a bit around to be harmonize with the NCBI metadata

metagenome_img =
  read_tsv("data/metadata_img.tsv",
          name_repair = "universal") %>%
  left_join(read_csv("data/IMG.csv",
                     name_repair = "universal") %>%
              select(IMG.Genome.ID,
                     site)) %>%
  mutate(n_genome = n(),
         .by = site) %>%
  select(accession = Bin.ID,
         n_genes = Gene.Count,
         bioproject = Genome.Name,
         n_bases = Total.Number.of.Bases,
         site,
         n_genome,
         GTDB.Taxonomy.Lineage) %>%
  separate_wider_delim(GTDB.Taxonomy.Lineage, delim = "; ",
                       names = c("superkingdom",
                                 "phylum",
                                 "class",
                                 "order",
                                 "family",
                                 "genus",
                                 "species"),
                                 cols_remove = T, too_few = "align_start")

# finally, we join IMG and NCBI metadata

metadata_metagenome =
  bind_rows(metagenome_ncbi, metagenome_img) %>%
  mutate(site = as.factor(site) %>%
           fct_lump_min(50, other_level = "Multiple Sites"))

count(metadata_metagenome, site)

# load tabular results from hmmsearch

hmm_raw_results =
  list.files("data/results_hmmer",
             pattern = ".csv", full.names = T) %>%
  map(~ read_csv(.x,
                 id = "accession",
                 col_names = c("target", "query", "evalue", "score", "bias",
                                "dom_evalue", "dom_score", "dom_bias"))) %>%
  list_rbind()

# add metadata on genomes and genes

results_hmmer =
  hmm_raw_results %>%
  mutate(accession = str_remove(accession, "data/results_hmmer/") %>%
           str_remove(".csv"),
         Locustag = str_remove(query, "-i1")) %>%
  left_join(metadata_metagenome) %>%
  left_join(dcgs)

# write this down cz the operation takes time

write_csv(results_hmmer, "results_hmmer.csv")
results_hmmer = read_csv("results_hmmer.csv")

# filtering likely spurious hits

res_filtered =
  results_hmmer %>%
  filter(
    dom_score > 0,
    score > 0,
    score > 2 * bias,
    dom_score > 2 * dom_bias,
    score < 2 * dom_score,
    dom_evalue < 1/n_genes)

res_filtered %>%
  count(Gene, site) %>% View()

# plot taxonomical distribution across sites

res_filtered %>%
  mutate(order2 =
           as.factor(order) %>%
           fct_na_value_to_level("Other") %>%
           fct_lump_n(5, other_level = "Other") %>%
           factor(levels = c("Other", "Anaerolineales",
                             "Halobacteriales", "Pseudomonadales",
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
  facet_grid(Operon ~ site, scales = "free",
             space = "free_y") +
  scale_fill_manual(values = c("grey95", "grey80", "grey65", "grey50", "#FFC107")) +
  scale_x_log10()


ggsave("metagenomes.svg", width = 24, height = 12)