#This script will investigate the occurence of bacteria in three datasets using three different appraoches and determine the shared bacteria that may lead to a reduced prognosis
print('Loading packages and printing package versions...')
library(phyloseq); print('phyloseq'); packageVersion('phyloseq')
library(ggplot2); print('ggplot2'); packageVersion('ggplot2')
library(plyr); print('plyr'); packageVersion('plyr')
library(grid); print('grid'); packageVersion('grid')
library(scales); print('scales'); packageVersion('scales')
library(RColorBrewer); print('RColorBrewer'); packageVersion('RColorBrewer')
library(ComplexHeatmap); print('ComplexHeatmap'); packageVersion('ComplexHeatmap')
library(circlize); print('circlize'); packageVersion('circlize')
library(ggrepel); print('ggrepel'); packageVersion('ggrepel')
library(mclust); print('mclust'); packageVersion('mclust')
library(reshape2); print('reshape2'); packageVersion('reshape2')
library(gridExtra); print('gridExtra'); packageVersion('gridExtra')
library(survminer); print('survminer'); packageVersion('survminer')
library(survival); print('survival'); packageVersion('survival')
library(ggforce); print('ggforce'); packageVersion('ggforce')
library(tidyverse); print('tidyverse'); packageVersion('tidyverse')
library(datapasta); print('datapasta'); packageVersion('datapasta')
library(ggpubr); print('ggpubr'); packageVersion('ggpubr')
library(cowplot); print('cowplot'); packageVersion('cowplot')
library(ggbeeswarm); print('ggbeeswarm'); packageVersion('ggbeeswarm')
library(corrr); print('corrr'); packageVersion('corrr')
library(datapasta); print('datapasta'); packageVersion('datapasta')
library(lubridate); print('lubridate'); packageVersion('lubridate')
library(kableExtra); print('kableExtra'); packageVersion('kableExtra')
library(readxl); print('readxl'); packageVersion('readxl')
library(vegan); print('vegan'); packageVersion('vegan')
library(svglite); print('svglite'); packageVersion('svglite')
library(lubridate)
library(magrittr)

setwd('~/Desktop/bacteria_repo/cross_dataset_survival/scripts/')

#set random seed
set.seed(98267459)
options(scipen=999)

#define pallette for survival data 
survival_palette <- c('firebrick4', 'deepskyblue4')

#define a function for km plots
make_km <- function(km_fit, input_data, xlimval = NULL) {
  
  if  (is.null(xlimval)) {
    xlimval = max(km_fit$time)
  }
  
  km_plot <- ggsurvplot(km_fit,
                             size=2,
                             pval.size=8,
                             data=input_data,
                             pval=TRUE,
                             legend.title='ABBS', legend.labs=c('Pos', 'Neg'),# A = not clustered (chartreuse),     B = clustered (darkorchid)
                             palette = survival_palette,
                             ggtheme = theme_pubclean(base_size=17),
                             tables.theme=theme_survminer(), break.time.by = 500,
                             conf.int = FALSE,
                             risk.table = TRUE,
                             cumevents=TRUE, censor=FALSE, #censor.size=0, #changed from 8 17th Feb 2022 following reviewer request
                             ylab = 'Progression Free Survival',
                             xlab='Days To Event',
                        xlim = c(0, xlimval))
  
 return(km_plot)
}
  

#define interesting genera to investigate
interesting_genera = c('Ezakiella', 'Sporobacterium', 'Fenollaria', 'Peptoniphilus', 'Porphyromonas', 'Anaerococcus', 'Fusobacterium')

#Load in community data frames for ICGC, 16s, RNAseq
icgc_community <- read_tsv(file='../data/ICGC_Tissue_Community.tsv', col_names = TRUE)
s16_community <- read_tsv(file='../data/16s_Community.tsv', col_names = TRUE)
rnaseq_community <- read_tsv(file='../data/rnaseq_Community.tsv', col_names=TRUE)

#Load in survival data for ICGC, 16s, RNAseq
icgc_survival <- read_tsv(file='../data/ICGC_Progression_Data.tsv', col_names = TRUE)
s16_survival <-read_tsv(file='../data/16s_Progression_Data.tsv', col_names=TRUE)

#rnaseq_survival <- read_tsv(file='../data/rnaseq_Progression_Data.tsv', col_names=TRUE)

#minus 1 from survival binary identifiers to match other data frames
s16_survival <- s16_survival %>% filter(progression_days >= 0) %>%
  mutate(progression = progression -1,
         skeletal_metastases = case_when(
           skeletal_metastases == 'yes' ~ 1,
           skeletal_metastases == 'no' ~ 0))

#calculate time to metastases in 16s
s16_survival <- s16_survival %>% mutate(enrollment_date = ymd(enrollment_date),
                                        followUp_date = ymd(followUp_date),
                                        mets_date = ymd(mets_date))

s16_survival <- s16_survival %>% mutate(
  mets_days = case_when(
    skeletal_metastases == 1 ~ (mets_date-enrollment_date),
    skeletal_metastases == 0 ~ (followUp_date-enrollment_date))
)

s16_survival$mets_days <- as.numeric(s16_survival$mets_days)

#mutate progression date to be the min of mets date and progression date -- because mets or progression is still progression
s16_survival <- s16_survival %>% mutate(
  progression_days = pmin(progression_days, mets_days)
  )



#### Start analysing survival with defined interesting Taxa - clean and format input data ####

#set min number of genera to investigate
min_num_genera <- 1

# Clean and format ICGC Data
#work out which patients had interesting genera
icgc_filt <- icgc_community %>% filter(Taxa %in% interesting_genera)
icgc_filt <- icgc_filt %>% column_to_rownames('Taxa')
icgc_filt <- as.data.frame(t(icgc_filt))
#isolate counts of interesting genera
icgc_filt$rowsum <- rowSums(icgc_filt) 
icgc_filt <- icgc_filt %>% rownames_to_column('SampleId')
#Mark interesting genera
icgc_filt <- icgc_filt %>% mutate(interesting_genera = case_when(
  rowsum >= min_num_genera ~ 'Interesting_Genera',
  TRUE ~ 'No_Interesting_Genera'
))

# find mets samples to remove
icgc_sample_data <- read_tsv(file='../data/ICGC_samples_w_batch_info.tsv')
icgc_sample_filt <- icgc_sample_data %>% select(Sanger_ID, ICGC_sample_ID, sample_type) %>%
  filter(sample_type != 'normal - blood')
mets_samples <- icgc_sample_filt %>% filter(sample_type %in% c('metastatic tissue (autopsy)',
                                                               'metastatic tissue (at prostatectomy)'))
# filter out mets samples -- keep only primary tumours - otherwise results would be biased for patients with multiple samples (mets)
icgc_filt <- icgc_filt %>% filter(!SampleId %in% mets_samples$Sanger_ID)




#Load in ID conversion file
icgc_id <- read_tsv(file='../data/ICGC_ID_Converstion.tsv', col_names=TRUE)
#remove mets samples from conversion
icgc_id <- icgc_id %>% filter(!SampleId %in% mets_samples$Sanger_ID)

icgc_filt <- merge(icgc_filt, icgc_id, by='SampleId', all.x=TRUE, all.y=FALSE)
icgc_filt <- icgc_filt %>% select(ICGC_ID, interesting_genera)
# merge in interesting_genera into survival data
icgc_survival <- merge(icgc_filt, icgc_survival, by='ICGC_ID', all=TRUE)
# filter out NA survival
icgc_survival <- icgc_survival %>% filter(!is.na(event_type))



# Clean and format RNA seq data

##presence absence and filter for only interesting taxa
#rnaseq_community[,2:ncol(rnaseq_community)] <- decostand(rnaseq_community[,2:ncol(rnaseq_community)], method='pa')
#rnaseq_community <- rnaseq_community %>% filter(name %in% interesting_genera)
#rnaseq_community <- rnaseq_community %>% column_to_rownames('name')
#rnaseq_community <- as.data.frame(t(rnaseq_community))
#rnaseq_community$rowsum <- rowSums(rnaseq_community)
#rnaseq_community <- rnaseq_community %>% rownames_to_column('name')
# 
# rnaseq_community <- rnaseq_community %>% mutate(interesting_genera = case_when(
#   rowsum >= min_num_genera  ~ 'Interesting_Genera',
#   TRUE ~ 'No_Interesting_Genera'
# ))
# 
# #merge in to survival
# rnaseq_community$name  <- gsub('P', '', rnaseq_community$name)
# rnaseq_community <- rnaseq_community %>% select(name,interesting_genera)
# colnames(rnaseq_community) <- c('id', 'interesting_genera')
# rnaseq_survival <- merge(rnaseq_survival, rnaseq_community, by='id', all=TRUE)



#  Clean and format 16s data
s16_community <- s16_community %>% column_to_rownames('sample_id')
s16_community$ids <- NULL

s16_names <- data.frame(taxa_string = colnames(s16_community))
s16_names <- separate(s16_names, taxa_string, c('Taxa', NULL), sep=';', remove=TRUE, convert=FALSE, extra='drop', fill='right')
s16_names <- separate(s16_names, Taxa, c('Extra', 'Taxa'), sep='__', remove=TRUE, convert=FALSE, extra='drop', fill='right')
#save end bit to add on to names later
end_bit <- s16_names$Extra[158:162]
s16_names <- s16_names %>% filter(!is.na(Taxa)) %>% select(Taxa)
# manually clean up extra bits with weird names
s16_names$Taxa <- gsub('Prevotella 6', 'Prevotella', s16_names$Taxa)
s16_names$Taxa <- gsub('Treponema 2', 'Treponema', s16_names$Taxa)
s16_names$Taxa <- gsub('Escherichia-Shigella', 'Escherichia', s16_names$Taxa)

reformatted_16s_names <- c(s16_names$Taxa, end_bit)

#manual check reformatted names match -- should work for majority of genera
for (i in seq(1,ncol(s16_community),1)) {
  print(paste0( colnames(s16_community)[i], ' --- ', reformatted_16s_names[i]  ))
}

colnames(s16_community) <- reformatted_16s_names
#filter out taxa < 5% abundance and make binary presence absence
s16_community[,1:157][s16_community[,1:157]  < 5 ] <- 0
s16_community[,1:157] <- decostand(s16_community[,1:157], method='pa')
#select only colnames that are in interesting_genera
s16_community <- s16_community %>% rownames_to_column('sample_ID')
s16_community <- s16_community %>% select(c(sample_ID, intersect(colnames(s16_community), interesting_genera), end_bit))

#count number of interesting genera in each sample -- excluding other factors
s16_community$rowsum <- rowSums(s16_community[,2:(ncol(s16_community) -5) ])

s16_selection <- s16_community %>% mutate(interesting_genera = case_when(
  rowsum >= min_num_genera  ~ 'Interesting_Genera',
  TRUE ~ 'No_Interesting_Genera'
)) %>% select(sample_ID, interesting_genera)

#merge into survival data and filter appropriately
s16_survival <- merge(s16_survival, s16_selection, by='sample_ID', all=TRUE)
s16_survival <- s16_survival %>% filter(!is.na(progression))



#filter out the samples that were included as UTI controls and those with no clinical follow-up
s16_survival <- s16_survival %>% filter(progression_days > 0) %>% 
  filter(!sample_ID %in% c('M138.4', 'M86.1', 'M86.6'))



#### now do some km analysis for icgc_survival, rnaseq_survival and s16_survival ####
#model survival data
icgc_km_fit <- survfit(Surv(days_to_event ,event_type) ~ interesting_genera, data=icgc_survival)
#plot the km data
icgc_km_plot <- make_km(icgc_km_fit, icgc_survival, xlimval = sort(icgc_survival$days_to_event)[nrow(icgc_survival)-10]) #surv_summary(icgc_km_fit) # point at which n=4 ABBS+ = 2876   | ABBS- = 3352
icgc_km_plot



### load in updated clinical data as provided by rachel -- line updated 14th Aug 2020
rnaseq_survival <- data.frame(
                         stringsAsFactors = FALSE,
                                                   sample_id = c("P100_3",
                                                                 "P104_1","P104_4",
                                                                 "P109_5","P110_1",
                                                                 "P112_4","P113_1",
                                                                 "P114_5","P115_1",
                                                                 "P117_7","P123_5",
                                                                 "P124_2","P125_5",
                                                                 "P129_4","P131_2",
                                                                 "P133_4","P134_5",
                                                                 "P137_1","P138_6",
                                                                 "P140_2","P140_4",
                                                                 "P140_6","P143_1",
                                                                 "P143_2","P143_3",
                                                                 "P147_1","P148_2",
                                                                 "P26_2","P33_6",
                                                                 "P45_2","P48_4",
                                                                 "P62_1","P75_2","P82_7",
                                                                 "P86_3","P88_7",
                                                                 "P97_7","P98_1",
                                                                 "P99_2","P99_4"),
                                               days_to_event = c(2116L,2053L,
                                                                 2053L,1636L,
                                                                 2024L,1946L,2010L,
                                                                 1953L,665L,821L,
                                                                 1069L,1906L,519L,
                                                                 1871L,1843L,1618L,
                                                                 870L,1793L,1787L,
                                                                 1779L,1779L,1779L,
                                                                 1758L,1758L,1758L,
                                                                 1737L,1731L,262L,
                                                                 2535L,2493L,2480L,
                                                                 1407L,2256L,2165L,
                                                                 1648L,2221L,878L,
                                                                 2144L,2137L,2137L),
                                                  event_type = c(0L,0L,0L,
                                                                 1L,0L,1L,0L,1L,
                                                                 0L,1L,1L,0L,1L,
                                                                 0L,0L,1L,1L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,1L,
                                                                 0L,0L,0L,1L,1L,
                                                                 0L,1L,0L,1L,0L,
                                                                 0L,0L),
                        category_at_initial_urine_collection = c("I","I",
                                                                 "H","I","CB","I",
                                                                 "CB","H","CB","H",
                                                                 "H","H","I","H",
                                                                 "CB","H","I",
                                                                 "CB","CB","CB","CB",
                                                                 "CB","CB","CB",
                                                                 "CB","CB","H","CB",
                                                                 "CB","CB","I",
                                                                 "I","H","I","H",
                                                                 "H","I","I","H",
                                                                 "I"),
                     subcategory_at_initial_urine_collection = c("I","I",
                                                                 "HL","I","CB<1","I",
                                                                 "CBN","Hh",
                                                                 "CB<1","Hh","Hh","Hh",
                                                                 "I","HL","CB<1",
                                                                 "Hh","IH","CBN",
                                                                 "CB<1","CBN","CB<1",
                                                                 "CBN","CBN",
                                                                 "CB<1","CBN","CBN",
                                                                 "Hh","CBN","CBN",
                                                                 "CBN","I","I","HL",
                                                                 "I","Hh","Hh",
                                                                 "I","I","Hh","I"),
                                                        MDS1 = c(-0.613935794,
                                                                 -0.897081171,
                                                                 -1.088960908,
                                                                 0.287979616,0.667874573,
                                                                 0.365655686,0.019834977,
                                                                 -1.115587403,
                                                                 0.578313482,0.377639909,
                                                                 0.719965531,
                                                                 -0.759494243,-0.26101531,
                                                                 -0.082201755,
                                                                 0.551544341,0.466636606,
                                                                 0.379863061,
                                                                 0.268458015,0.377953495,
                                                                 -0.213329126,
                                                                 0.533168414,0.332058631,
                                                                 0.160642049,
                                                                 0.623366707,0.540049411,
                                                                 -1.01505389,
                                                                 -0.745882232,-0.30254766,
                                                                 -0.807282442,
                                                                 -0.22873891,0.18540428,
                                                                 -0.322742219,
                                                                 0.679668286,-0.381623413,
                                                                 0.596919072,
                                                                 0.565893118,0.225128976,
                                                                 -0.753363681,
                                                                 -0.576869861,0.661691784),
                                                        MDS2 = c(0.803296633,
                                                                 -0.732098322,
                                                                 0.081733134,0.127440513,
                                                                 0.052577932,
                                                                 0.35441864,0.011302458,
                                                                 0.401244978,
                                                                 0.295582447,0.484772383,
                                                                 -0.098046263,
                                                                 0.555472881,-0.566415923,
                                                                 -0.253284886,
                                                                 0.188661125,-0.205725277,
                                                                 -0.581940327,
                                                                 0.52578929,0.2162203,
                                                                 0.287835081,
                                                                 0.149771671,-0.089182162,
                                                                 0.096833956,
                                                                 -0.012262512,0.083479043,
                                                                 0.57545525,
                                                                 -1.030093334,-0.220864052,
                                                                 -0.267601539,
                                                                 0.348640044,-0.52511186,
                                                                 0.105600324,
                                                                 -0.477547333,-0.045366729,
                                                                 0.57858378,
                                                                 -0.205218904,-0.34565192,
                                                                 -0.124898222,
                                                                 -0.205839342,-0.337562956),
                                                     Cluster = c(2L,2L,2L,
                                                                 1L,1L,1L,1L,2L,
                                                                 1L,1L,1L,2L,2L,
                                                                 1L,1L,1L,1L,1L,
                                                                 1L,2L,1L,1L,1L,
                                                                 1L,1L,2L,2L,2L,
                                                                 2L,2L,1L,2L,1L,
                                                                 2L,1L,1L,1L,2L,
                                                                 2L,1L),
                                               cluster_label = c("A","A",
                                                                 "A","B","B","B",
                                                                 "B","A","B","B",
                                                                 "B","A","A","B",
                                                                 "B","B","B","B",
                                                                 "B","A","B","B",
                                                                 "B","B","B","A",
                                                                 "A","A","A","A",
                                                                 "B","A","B","A",
                                                                 "B","B","B","A",
                                                                 "A","B"),
                                                   Ezakiella = c(0L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,164L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,158L),
                                               Peptoniphilus = c(0L,0L,0L,
                                                                 0L,360L,0L,0L,0L,
                                                                 0L,187L,334L,0L,
                                                                 0L,0L,0L,103L,
                                                                 1536L,0L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,215L,5710L,0L,
                                                                 0L,289L,139L,0L,
                                                                 0L,1420L),
                                               Porphyromonas = c(0L,0L,0L,
                                                                 0L,16147L,0L,0L,
                                                                 0L,0L,0L,1473L,
                                                                 0L,0L,0L,113L,
                                                                 1672L,0L,0L,0L,0L,
                                                                 0L,0L,0L,155L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,76946L,0L,
                                                                 0L,10632L,0L,0L,
                                                                 0L,30872L),
                                                Anaerococcus = c(0L,0L,0L,
                                                                 0L,288L,0L,0L,0L,
                                                                 0L,0L,498L,0L,
                                                                 129L,0L,0L,137L,
                                                                 372L,0L,0L,0L,0L,
                                                                 0L,0L,0L,1055L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,455L,0L,
                                                                 0L,276L,153L,0L,
                                                                 0L,531L),
                                               Fusobacterium = c(0L,0L,0L,
                                                                 183L,119L,184L,0L,
                                                                 0L,0L,148L,0L,
                                                                 0L,0L,0L,271L,
                                                                 222L,0L,0L,0L,0L,
                                                                 0L,0L,192L,0L,
                                                                 125L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,139L,
                                                                 267L,2553L,0L,
                                                                 0L,0L,2119L),
                                               abbs_read_sum = c(0L,0L,0L,
                                                                 183L,16914L,184L,
                                                                 0L,0L,0L,335L,
                                                                 2305L,0L,129L,0L,
                                                                 384L,2134L,2072L,
                                                                 0L,0L,0L,0L,0L,
                                                                 192L,155L,1180L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 215L,83111L,139L,
                                                                 267L,13750L,292L,
                                                                 0L,0L,35100L),
                                                  abbs_group = c("Non-ABBS",
                                                                 "Non-ABBS",
                                                                 "Non-ABBS","ABBS","ABBS",
                                                                 "ABBS","Non-ABBS",
                                                                 "Non-ABBS",
                                                                 "Non-ABBS","ABBS","ABBS",
                                                                 "Non-ABBS","ABBS",
                                                                 "Non-ABBS","ABBS",
                                                                 "ABBS","ABBS",
                                                                 "Non-ABBS","Non-ABBS",
                                                                 "Non-ABBS",
                                                                 "Non-ABBS","Non-ABBS",
                                                                 "ABBS","ABBS","ABBS",
                                                                 "Non-ABBS",
                                                                 "Non-ABBS","Non-ABBS",
                                                                 "Non-ABBS","Non-ABBS",
                                                                 "Non-ABBS","ABBS",
                                                                 "ABBS","ABBS",
                                                                 "ABBS","ABBS","ABBS",
                                                                 "Non-ABBS",
                                                                 "Non-ABBS","ABBS"),
                                                strict_event = c(0L,0L,0L,
                                                                 1L,0L,1L,0L,1L,
                                                                 0L,1L,1L,0L,1L,
                                                                 0L,0L,1L,0L,0L,
                                                                 0L,0L,0L,0L,0L,
                                                                 0L,0L,0L,0L,1L,
                                                                 0L,0L,0L,1L,1L,
                                                                 0L,1L,0L,1L,0L,
                                                                 0L,0L),
                                           strict_event_date = c(NA,NA,NA,
                                                                 "22/09/2014",NA,
                                                                 "03/03/2014",NA,
                                                                 "27/12/2018",NA,
                                                                 "12/09/2016","16/04/2019",
                                                                 NA,"23/08/2015",
                                                                 NA,NA,"15/12/2018",
                                                                 "22/09/2019",NA,
                                                                 NA,NA,NA,NA,NA,
                                                                 NA,NA,NA,NA,
                                                                 "20/05/2016",NA,NA,NA,
                                                                 "13/11/2018",
                                                                 "18/02/2020",NA,
                                                                 "22/09/2017",NA,
                                                                 "13/11/2015",NA,NA,NA),
                                           strict_event_type = c(NA,NA,NA,
                                                                 "PSA_prog",NA,
                                                                 "PSA_prog",NA,"BCR",
                                                                 NA,"mets","mets",
                                                                 NA,"PSA_prog",NA,
                                                                 NA,"BCR","RP",NA,
                                                                 NA,NA,NA,NA,NA,
                                                                 NA,NA,NA,NA,
                                                                 "BCR",NA,NA,NA,"BCR",
                                                                 "BCR",NA,"BCR",
                                                                 NA,"PSA_prog",NA,
                                                                 NA,NA),
                                            recruitment_date = c(NA,NA,NA,
                                                                 "08/10/2013",NA,
                                                                 "22/10/2013",NA,
                                                                 "29/10/2013",NA,
                                                                 "12/11/2013","28/01/2014",
                                                                 NA,"11/02/2014",
                                                                 NA,NA,"29/04/2014",
                                                                 "06/05/2014",NA,
                                                                 NA,NA,NA,NA,NA,
                                                                 NA,NA,NA,NA,
                                                                 "24/04/2012",NA,NA,NA,
                                                                 "25/09/2012",
                                                                 "11/12/2012",NA,
                                                                 "19/03/2013",NA,
                                                                 "18/06/2013",NA,NA,NA)
                   )

rnaseq_survival$strict_event_date <- dmy(rnaseq_survival$strict_event_date)
rnaseq_survival$recruitment_date <- dmy(rnaseq_survival$recruitment_date)

rnaseq_survival$strict_days_to_event <- rnaseq_survival$strict_event_date - rnaseq_survival$recruitment_date

#use last followup date if no updated event
rnaseq_survival %<>% mutate(strict_days_to_event_complete = case_when(
  as.numeric(strict_days_to_event) > -1  ~ as.double(strict_days_to_event),
  TRUE ~ as.double(days_to_event)
))

#write.table(rnaseq_survival, file='../data/rnaseq_cancer_survival_RNA_Strict.tsv', col.names = TRUE, row.names = FALSE, quote=FALSE, sep='\t')


#model survival data
rnaseq_km_fit <- survfit(Surv(strict_days_to_event_complete, strict_event) ~ abbs_group, data=rnaseq_survival)
#plot the km data
rnaseq_km_plot <- make_km(rnaseq_km_fit, rnaseq_survival, xlimval=sort(rnaseq_survival$days_to_event)[nrow(rnaseq_survival)-10])
rnaseq_km_plot


#rnaseq only cancer samples 
rna_surv <- rnaseq_survival %>% filter(strict_event==1 | category_at_initial_urine_collection %in% c("A", "H", "I", "L"))
rna_km <-  survfit(Surv(strict_days_to_event_complete ,strict_event) ~ abbs_group, data=rna_surv)
rnaseq_km_p <- make_km(rna_km,rna_surv, xlimval = sort(rna_surv$days_to_event)[nrow(rna_surv)-10]) # point at which n=4 ABBS+ = 2165 | ABBS- = 2116
rnaseq_km_p


#model survival data
s16_km_fit <- survfit(Surv(progression_days , progression) ~ interesting_genera, data=s16_survival)
#plot the km data
s16_km_plot <- make_km(s16_km_fit,s16_survival)
s16_km_plot



#16s analysis with only evidence of cancer
s16_surv <- s16_survival %>% filter(progression==1 | category  %in% c("A", "H", "I", "L"))
s16_km <- survfit(Surv(progression_days , progression) ~ interesting_genera, data=s16_surv)
#plot the km data with only cancer samples
s16_km_p <- make_km(s16_km, s16_surv, xlimval = sort(s16_surv$progression_days)[nrow(s16_surv)-10]) # point at which n=4 ABBS+ = 945   | ABBS- = 839
s16_km_p



#16s mets
#model survival data
s16_km_fit_mets <- survfit(Surv(mets_days , skeletal_metastases) ~ interesting_genera, data=s16_survival)
#plot the km data
s16_km_plot_mets <- make_km(s16_km_fit_mets, s16_survival)
s16_km_plot_mets

#16s analysis with only evidence of cancer
s16_surv_mets <- s16_survival %>% filter(progression==1 | category  %in% c("A", "H", "I", "L"))
s16_km_mets <- survfit(Surv(mets_days , skeletal_metastases) ~ interesting_genera, data=s16_surv_mets)
#plot the km data with only cancer samples
s16_km_p_mets <- make_km(s16_km_mets, s16_surv_mets)
s16_km_p_mets




#Produce a nice object from the survival plot
icgc_ggplot <- ggdraw() +
  draw_plot(icgc_km_plot$plot, x=0,y=0.3,width=1,height=0.65) +
  draw_plot(icgc_km_plot$table, x=0, y=0.15, width=1, height=0.15 ) +
  draw_plot(icgc_km_plot$cumevents, x=0, y=0, width=1, height=0.15)
rnaseq_ggplot <- ggdraw() +
  draw_plot(rnaseq_km_plot$plot, x=0,y=0.3,width=1,height=0.65) +
  draw_plot(rnaseq_km_plot$table, x=0, y=0.15, width=1, height=0.15 ) +
  draw_plot(rnaseq_km_plot$cumevents, x=0, y=0, width=1, height=0.15)
s16_ggplot <- ggdraw() +
  draw_plot(s16_km_plot$plot, x=0, y=0.3,width=1,height=0.65) +
  draw_plot(s16_km_plot$table, x=0, y=0.15, width=1, height=0.15 ) +
  draw_plot(s16_km_plot$cumevents, x=0, y=0, width=1, height=0.15)

s16_mets_ggplot <- ggdraw() +
  draw_plot(s16_km_plot_mets$plot, x=0, y=0.3,width=1,height=0.65) +
  draw_plot(s16_km_plot_mets$table, x=0, y=0.15, width=1, height=0.15 ) +
  draw_plot(s16_km_plot_mets$cumevents, x=0, y=0, width=1, height=0.15)
 
 
#just cancer samples
s16_km_ggp <- ggdraw() +
  draw_plot(s16_km_p$plot, x=0, y=0.3,width=1,height=0.65) +
  draw_plot(s16_km_p$table, x=0, y=0.15, width=1, height=0.15 ) +
  draw_plot(s16_km_p$cumevents, x=0, y=0, width=1, height=0.15)

s16_km_ggp_mets <- ggdraw() +
  draw_plot(s16_km_p_mets$plot, x=0, y=0.3,width=1,height=0.65) +
  draw_plot(s16_km_p_mets$table, x=0, y=0.15, width=1, height=0.15 ) +
  draw_plot(s16_km_p_mets$cumevents, x=0, y=0, width=1, height=0.15)

rnaseq_km_ggp <- ggdraw() +
  draw_plot(rnaseq_km_p$plot, x=0, y=0.3,width=1,height=0.65) +
  draw_plot(rnaseq_km_p$table, x=0, y=0.15, width=1, height=0.15 ) +
  draw_plot(rnaseq_km_p$cumevents, x=0, y=0, width=1, height=0.15)

#print plots
icgc_ggplot
rnaseq_ggplot
s16_ggplot

#just cancer free plots
s16_km_ggp
rnaseq_km_ggp


km_plots <- cowplot::plot_grid(s16_ggplot, rnaseq_ggplot, icgc_ggplot, labels=c('a) 16s', 'b) RNA Sequencing', 'c) ICGC'), 
                               nrow=1, scale=0.9, label_size = 17)
km_plots
cancer_km_plots <- cowplot::plot_grid(s16_km_ggp, rnaseq_km_ggp, icgc_ggplot, labels=c('a) 16s', 'b) RNA Sequencing', 'c) ICGC'),
                                      nrow=1, scale=0.9, label_size=17)
cancer_km_plots


km_cancer='../plots/km_just_cancer.pdf'
cancer_svg = gsub('.pdf', '.svg', km_cancer)
ggsave(cancer_km_plots, file=km_cancer, width=25, height=10)
ggsave(cancer_km_plots, file=cancer_svg, width=25, height=10)
opencommand=paste0('open ', km_cancer)
system(opencommand)

kmfilename='../plots/km_all_interesting_taxa_1_minimum.pdf'
#ggsave(km_plots, file=kmfilename, width=20, height=10)
opencommand = paste0('open ', kmfilename)
system(opencommand)


# load previous 16s plots

combined_plot <- readRDS(file='../data/16s_plots.RDS')



new_cancer_km_plots <- cowplot::plot_grid(s16_km_ggp, rnaseq_km_ggp, icgc_ggplot, labels=c('d', 'e', 'f'),
                                      nrow=1, scale=0.9, label_size=30)


figure2 <- plot_grid(combined_plot, new_cancer_km_plots, labels=NULL, nrow=2)
#figure2

ggsave(figure2, filename = '../plots/new_fig2_strict_rna.pdf', height=20, width=20)
ggsave(figure2, filename = '../plots/new_fig2_strict_rna.svg', height=20, width=20)


#save survival dataframes for coxph later (only cancer samples)
# system('mkdir -pv ../data/cancer_samples')
# write.table(rna_surv, file='../data/cancer_samples/rnaseq_cancer_survival.tsv', col.names = TRUE, row.names = FALSE, sep='\t')
# write.table(s16_surv, file='../data/cancer_samples/16s_cancer_survival.tsv', col.names=TRUE, row.names=FALSE, sep='\t')
# write.table(icgc_survival, file='../data/cancer_samples/icgc_survival.tsv', col.names = TRUE, row.names = FALSE, sep='\t')



icgc_without_event <- icgc_survival %>% filter(event_type == 0)
s16_surv_no_event <- s16_surv %>% filter(progression==0)
rna_surv_no_event <- rna_surv %>% filter(event_type ==0)


summary(s16_surv_no_event$progression_days)


# lok int total number of events in each group
icgc_bac <- icgc_survival %>% filter(interesting_genera == 'Interesting_Genera')
icgc_no_bac <- icgc_survival %>% filter(interesting_genera == 'No_Interesting_Genera')

sum(icgc_survival$event_type)
sum(icgc_bac$event_type)
sum(icgc_no_bac$event_type)


s16_bac <- s16_surv %>% filter(interesting_genera == 'Interesting_Genera')
s16_no_bac <- s16_surv %>% filter(interesting_genera == 'No_Interesting_Genera')

sum(s16_bac$progression)
sum(s16_no_bac$progression)


rna_bac <- rna_surv %>% filter(abbs_group == 'ABBS')
rna_no_bac <- rna_surv %>% filter(abbs_group=='Non-ABBS')

sum(rna_bac$event_type)
sum(rna_no_bac$event_type)






