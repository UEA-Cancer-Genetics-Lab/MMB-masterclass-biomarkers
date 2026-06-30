library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(ggsurvfit)
library(survival)
library(Boruta)
library(DataExplorer)

s16_clin_data <- read_delim("Input/16s_cancer_survival.tsv")
s16_clin_data %>% select(sample_ID, category, progression_days,progression)
s16_clin_data <- s16_clin_data %>% select(sample_ID, category, progression_days,progression)
s16_clin_data$category[s16_clin_data$category=="S"] <- "L"
s16_clin_data$category <- factor(s16_clin_data$category, levels=c("L","I","H","A"))

survfit2(Surv(progression_days, progression)~1, data=s16_clin_data) %>% ggsurvfit()

s16_community <- read_delim("Input/16s_Community.tsv")[,1:158]

s16_names <- data.frame(taxa_string = colnames(s16_community))
s16_names <- separate(s16_names, taxa_string, c('Taxa', NULL), sep=';', remove=TRUE, convert=FALSE, extra='drop', fill='right')
s16_names <- separate(s16_names, Taxa, c('Extra', 'Taxa'), sep='__', remove=TRUE, convert=FALSE, extra='drop', fill='right')

s16_names <- s16_names %>% filter(!is.na(Taxa)) %>% select(Taxa)
# manually clean up extra bits with weird names
s16_names$Taxa <- gsub('Prevotella 6', 'Prevotella', s16_names$Taxa)
s16_names$Taxa <- gsub('Treponema 2', 'Treponema', s16_names$Taxa)
s16_names$Taxa <- gsub('Escherichia-Shigella', 'Escherichia', s16_names$Taxa)

s16_names$Taxa[duplicated(s16_names$Taxa)] <- paste0(s16_names$Taxa[duplicated(s16_names$Taxa)],".2")

colnames(s16_community) <- c("sample_ID",s16_names$Taxa)

# same samples in both tables
common_ids <- intersect(s16_community$sample_ID, s16_clin_data$sample_ID)

s16_community <- s16_community %>% filter(sample_ID %in% common_ids)
s16_clin_data <- s16_clin_data %>% filter(sample_ID %in% common_ids)

save(s16_community, s16_clin_data, file = "16s_PCa_data.RData")

load("16s_PCa_data.RData")

# EDA
create_report(s16_community, output_file = "community_report.html")
create_report(s16_clin_data, output_file = "clin_data_report.html")

# Remove those values that are less than 5% and convert to presence/absence
s16_community <- s16_community %>% mutate_if(is.numeric, ~1 * (. > 5))

# Only select genera with more than 2 hits
s16_community <- s16_community %>% select_if(function(col) is.character(col) || (is.numeric(col) && sum(col) >2))

# Merge taxa and survival
s16_merge <-  s16_clin_data %>% left_join(s16_community)

# Peptoniphilus 
survfit2(Surv(progression_days, progression)~Peptoniphilus, data=s16_merge) %>% ggsurvfit() + add_pvalue()

## Stretch make production quality figure

# ABBS
abbs_genera <- c("Ezakiella","Peptoniphilus","Porphyromonas","Anaerococcus","Fusobacterium")
s16_merge$abbs <- rowSums(s16_merge[,colnames(s16_merge) %in% abbs_genera]) > 0

survfit2(Surv(progression_days, progression)~abbs, data=s16_merge) %>% ggsurvfit()

survdiff(Surv(progression_days, progression)~abbs, data=s16_merge)$pvalue

coxph(Surv(progression_days, progression)~abbs+(category=="A"), data=s16_merge)
coxph(Surv(progression_days, progression)~abbs+I(category %in% c("H","A")), data=s16_merge)

coxph(Surv(progression_days, progression)~abbs+I(category %in% c("H","A")), data=s16_merge) %>% tidy(exponentiate=TRUE, conf.int=TRUE)      


# Boruta
set.seed(1000)
boruta_res <-  Boruta(Surv(s16_merge$progression_days, s16_merge$progression)~.,data=s16_merge[,-1:-4] )
boruta_res
plot(boruta_res)

#RNAseq data
rnaseq_clin_data <- read_delim("Input/rnaseq_cancer_survival.tsv") %>%
  rename(sample_ID=id, category = category_at_initial_urine_collection, progression = event_type, progression_days = days_to_event)
rnaseq_clin_data <- rnaseq_clin_data %>% select(sample_ID, category, progression_days,progression)
rnaseq_clin_data$category <- factor(rnaseq_clin_data$category, levels=c("I","H"))
rnaseq_clin_data$sample_ID <- paste0("P",rnaseq_clin_data$sample_ID)

survfit2(Surv(progression_days, progression)~1, data=rnaseq_clin_data) %>% ggsurvfit()

rnaseq_community <- read_delim("Input/rnaseq_Community.tsv")

# same samples in both tables
common_ids <- intersect(rnaseq_community$sample_ID, rnaseq_clin_data$sample_ID)

rnaseq_community <- rnaseq_community %>% filter(sample_ID %in% common_ids)
rnaseq_clin_data <- rnaseq_clin_data %>% filter(sample_ID %in% common_ids)

save(rnaseq_community, rnaseq_clin_data, file = "rnaseq_PCa_data.RData")