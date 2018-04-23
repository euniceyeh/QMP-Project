# longitudinal time-series plots across all CD vs Control patients in HMP2 for Prevotella and E. coli (genus-level and species-level)

library(tidyverse)

# import RMP for all CD patients

RMP_CD <- read.table(file = "../data/derived/taxonomic_profiles_metagenomics_CD_subset_week_num_transposed.tsv", sep = '\t', header = TRUE) %>%
  select(X, wk = week_num, Prevotella = ends_with("g__Prevotella"), Escherichia = ends_with("g__Escherichia"),     # genus level
         P_copri = ends_with("s__Prevotella_copri"), E_coli = ends_with("s__Escherichia_coli")) %>% # species level
  mutate(Diagnosis = rep("CD")) # all are CD

# then RMP for all controls

RMP_con <- read.table(file = "../data/derived/taxonomic_profiles_metagenomics_nonIBD_subset_week_num_transposed.tsv", sep = '\t', header = TRUE) %>%
  select(X, wk = week_num, Prevotella = ends_with("g__Prevotella"), Escherichia = ends_with("g__Escherichia"),     # genus level
         P_copri = ends_with("s__Prevotella_copri"), E_coli = ends_with("s__Escherichia_coli")) %>% # species level
  mutate(Diagnosis = rep("Control")) # all are controls

# combine the two RMP datasets into one

RMP <- bind_rows(RMP_CD, RMP_con)

# check which weeks are relevant
table(RMP$Diagnosis, RMP$wk)

# calculate averages

RMP_avgs <- RMP %>% group_by(wk, Diagnosis) %>%
  summarise(Prevotella_avg = mean(Prevotella), P_copri_avg = mean(P_copri),
            Escherichia_avg = mean(Escherichia), E_coli_avg = mean(E_coli))

# Prevotella
# genus level
RMP_avgs %>% ggplot(aes(wk, Prevotella_avg, col = Diagnosis)) +
  geom_line(size=1) +
  # scale_y_continuous(trans = "log10") +
  ylab("Average Relative Genus Abundance (RMP)") +
  xlab("Week") + ggtitle("Prevotella") +
  theme_bw()
# species level
RMP_avgs %>% ggplot(aes(wk, P_copri_avg, col = Diagnosis)) +
  geom_line(size=1) +
  # scale_y_continuous(trans = "log10") +
  ylab("Average Relative Species Abundance (RMP)") +
  xlab("Week") + ggtitle("Prevotella copri") +
  theme_bw()

# E. coli
# genus level
RMP_avgs %>% ggplot(aes(wk, Escherichia_avg, col = Diagnosis)) +
  geom_line(size=1) +
  # scale_y_continuous(trans = "log10") +
  ylab("Average Relative Genus Abundance (RMP)") +
  xlab("Week") + ggtitle("Escherichia") +
  theme_bw()
# species level
RMP_avgs %>% ggplot(aes(wk, E_coli_avg, col = Diagnosis)) +
  geom_line(size=1) +
  # scale_y_continuous(trans = "log10") +
  ylab("Average Relative Species Abundance (RMP)") +
  xlab("Week") + ggtitle("Escherichia coli") +
  theme_bw()

# save images by 900 w x 500 h