install.packages("ggpval")
library(ggpval)
install.packages("ggpubr")
library(ggpubr)
library(ggplot2)
install.packages("rstatix")
library(rstatix)
install.packages("ggprism")
library(ggprism)
#making boxplots for SRR_117
##variants parsing


##single-end samples

filenames_single <- list.files("~/Desktop/tel_variants/single/computel_117_untrimmed_single_no_ref", pattern="SRR*", full.names=FALSE)
variants_single <- list.files("~/Desktop/tel_variants/single/computel_117_untrimmed_single_no_ref", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
print(filenames_single)
print(variants_single)
single <- lapply(variants_single, read.delim2)
single <- lapply(single, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(single)
single_subset <- lapply(single, function(x) subset(x, pattern != "TTAGGG"))
print(single_subset)
column_names <- c("pattern", "abs.num", "variants.percent")
single_subset <- lapply(single_subset, setNames, column_names)
names(single_subset) <- filenames_single
print(single_subset)

##paired samples
##computel_117_trimmed_no_ref
filenames_117_trimmed <- list.files("~/Desktop/tel_variants/paired/computel_117_trimmed_no_ref/output_1", pattern="SRR*", full.names=FALSE)
print(filenames_117_trimmed)
##output_1
variants_117_trimmed_1 <- list.files("~/Desktop/tel_variants/paired/computel_117_trimmed_no_ref/output_1", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
print(variants_117_trimmed_1)
trimmed_117_1 <- lapply(variants_117_trimmed_1, read.delim2)
trimmed_117_1 <- lapply(trimmed_117_1, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(trimmed_117_1)
names(trimmed_117_1) <- filenames_117_trimmed
##output_2
variants_117_trimmed_2 <- list.files("~/Desktop/tel_variants/paired/computel_117_trimmed_no_ref/output_2", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
trimmed_117_2 <- lapply(variants_117_trimmed_2, read.delim2)
trimmed_117_2 <- lapply(trimmed_117_2, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(variants_117_trimmed_2)
print(trimmed_117_2)
names(trimmed_117_2) <- filenames_117_trimmed

#combining output_1 and output_2
m <- match(names(trimmed_117_2), names(trimmed_117_1), nomatch = 0L)

print(m)
combined_117_trimmed_1_2 <- Map(rbind, trimmed_117_1[m], trimmed_117_2)
print(combined)

combined_117_trimmed_1_2 <- lapply(combined_117_trimmed_1_2, function(x) summarise_at(group_by(x, pattern), vars(abs.num), list(~floor(mean(.,na.rm=FALSE)))))
print(combined_117_trimmed_1_2)

percent_of_variants <- lapply(combined_117_trimmed_1_2, function(x) x[,sapply(x, is.numeric)]/colSums(x[,sapply(x, is.numeric)])*100)
print(percent_of_variants)
colname <- ("variants.percent")
percent_of_variants <- lapply(percent_of_variants, setNames, colname)
print(percent_of_variants)
combined_117_trimmed_1_2 <- Map(cbind, combined_117_trimmed_1_2, percent_of_variants) 
print(combined_117_trimmed_1_2)
combined_117_trimmed_1_2 <- lapply(combined_117_trimmed_1_2, function(x){
  x[order(x$abs.num, decreasing = T),]
})
print(combined_117_trimmed_1_2)
#combined_117_trimmed_1_2_top_5 <- lapply(combined_117_trimmed_1_2, function(x) x[c(2:6), ])
#combined_117_trimmed_1_2_subset <- lapply(combined_117_trimmed_1_2, function(x) subset(x, variants.percent > 1))
combined_117_trimmed_1_2_subset <- lapply(combined_117_trimmed_1_2, function(x) subset(x, pattern != "TTAGGG"))
print(combined_117_trimmed_1_2_subset)


##computel_117_untrimmed_no_ref
filenames_117_untrimmed <- list.files("~/Desktop/tel_variants/paired/computel_117_untrimmed_no_ref/output_1", pattern="SRR*", full.names=FALSE)
print(filenames_117_untrimmed)
##output_1
variants_117_untrimmed_1 <- list.files("~/Desktop/tel_variants/paired/computel_117_untrimmed_no_ref/output_1", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
print(variants_117_untrimmed_1)
untrimmed_117_1 <- lapply(variants_117_untrimmed_1, read.delim2)
untrimmed_117_1 <- lapply(untrimmed_117_1, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(untrimmed_117_1)
names(untrimmed_117_1) <- filenames_117_untrimmed
##output_2
variants_117_untrimmed_2 <- list.files("~/Desktop/tel_variants/paired/computel_117_untrimmed_no_ref/output_2", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
untrimmed_117_2 <- lapply(variants_117_untrimmed_2, read.delim2)
untrimmed_117_2 <- lapply(untrimmed_117_2, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(variants_117_untrimmed_2)
print(untrimmed_117_2)
names(untrimmed_117_2) <- filenames_117_untrimmed

#combining output_1 and output_2
q <- match(names(untrimmed_117_2), names(untrimmed_117_1), nomatch = 0L)
print(q)
combined_117_untrimmed_1_2 <- Map(rbind, untrimmed_117_1[q], untrimmed_117_2)
print(combined_117_untrimmed_1_2)

combined_117_untrimmed_1_2 <- lapply(combined_117_untrimmed_1_2, function(x) summarise_at(group_by(x, pattern), vars(abs.num), list(~floor(mean(.,na.rm=FALSE)))))
df <- lapply(combined_117_untrimmed_1_2, function(x) as.data.frame(x))
print(df)
print(combined_117_untrimmed_1_2)

percent_of_variants <- lapply(combined_117_untrimmed_1_2, function(x) x[,sapply(x, is.numeric)]/colSums(x[,sapply(x, is.numeric)])*100)
print(percent_of_variants)
colname <- ("variants.percent")
percent_of_variants <- lapply(percent_of_variants, setNames, colname)
print(percent_of_variants)
combined_117_untrimmed_1_2 <- Map(cbind, combined_117_untrimmed_1_2, percent_of_variants) 
print(combined_117_untrimmed_1_2)
combined_117_untrimmed_1_2 <- lapply(combined_117_untrimmed_1_2, function(x){
  x[order(x$abs.num, decreasing = T),]
})
print(combined_117_untrimmed_1_2)
#combined_117_untrimmed_1_2_top_5 <- lapply(combined_117_untrimmed_1_2, function(x) x[c(2:6), ])
#combined_117_untrimmed_1_2_subset <- lapply(combined_117_untrimmed_1_2, function(x) subset(x, variants.percent > 1))
combined_117_untrimmed_1_2_subset <- lapply(combined_117_untrimmed_1_2, function(x) subset(x, pattern != "TTAGGG"))
print(combined_117_untrimmed_1_2_subset)


#combining 117 everything
#variants_117 <- list(combined_117_untrimmed_1_2_top_5, combined_117_trimmed_1_2_top_5, single_top_5)
variants_117 <- list(combined_117_untrimmed_1_2_subset, combined_117_trimmed_1_2_subset, single_subset)
View(variants_117)
print(variants_117)
variants_117 <- do.call(c, variants_117)
variants_117 <- lapply(variants_117, function(x) x[(names(x) %in% c("pattern", "variants.percent"))])
View(variants_117)

tissue <- list_merge(variants_117[c(1, 19, 5, 20, 21, 22, 23, 24)])
cfDNA <- list_merge(variants_117[-c(1, 19, 5, 20, 21, 22, 23, 24)])
View(cfDNA)
SRR_117_cancer_healthy <- list(tissue, cfDNA)
names(SRR_117_cancer_healthy) <- c("tissue", "cfDNA")
print(SRR_117_cancer_healthy)

melt_117 <- melt(SRR_117_cancer_healthy, id = "pattern", value.name = "variants.percent")
View(melt_117)

#should I filter before or after grouping by?
#filtering before grouping by

melt_117_filtered <- subset(melt_117, variants.percent > 1)
View(melt_117_filtered)

#grouping by repeated variants to get a table with possible variants (>1%)
melt_117_filtered_variants <- melt_117_filtered %>% group_by(pattern, L1) %>%
  summarise(variants.percent=(mean(as.numeric(variants.percent))))
View(melt_117_filtered_variants)

#dealing with duplicates, this way you would just have variants with a percent in samples > 1%
melt_117_filtered_variants_no_dupl <- distinct(melt_117_filtered_variants, pattern)
View(melt_117_filtered_variants_no_dupl)

#calculating the average over all the samples, including ones with a percent < 1
melt_117_average <- melt_117[melt_117$pattern %in% melt_117_filtered_variants_no_dupl$pattern,] %>% group_by(pattern, L1) %>% summarise(variants.percent=(mean(as.numeric(variants.percent))))
View(melt_117_average)

#finding max in a table to put on y-axis
max_117 <- max(melt_117_new[melt_117_new$pattern %in% melt_117_filtered_variants_no_dupl$pattern,][, "variants.percent"])
print(max_117)

View(melt_117_new)
#a different approach with zeros
###the most common variants for cancer/healthy, tissue/cfDNA

#SRR_117

#there is a need to put zeros in variants percent if the certain variant doesn't exist in a sample
#if there is no variant from the melt_117_grouped for a sample, then set it to zero
#tissue samples
missing_values_tissue <- lapply(tissue, function(x) anti_join(melt_117_filtered_variants_no_dupl[,"pattern"], x, by = 'pattern'))
missing_values_tissue <- lapply(missing_values_tissue, cbind, variants.percent = c("0"))
print(missing_values_tissue)

tissue_new <- Map(rbind, tissue, missing_values_tissue)
View(tissue_new)
print(tissue_new)

#cfDNA samples
#creating functions that would find logical types in a list of dataframes and replace it with character
#logical types happen when the dataframe is empty
outer <- function(list_of_lists) {
  return (lapply(list_of_lists, inner))
}	

inner <- function(list) {
  return (data.frame(lapply(list, function(x) if(is.logical(x)) { 
    return(as.character(x))
  } else {  
    return(x)
  } )))
}	

cfDNA_character <- outer(cfDNA)

View(melt_117_grouped)
print(cfDNA_character)

missing_values_cfDNA <- lapply(cfDNA_character, function(x) anti_join(melt_117_filtered_variants_no_dupl[,"pattern"], x, by = 'pattern'))

missing_values_cfDNA <- lapply(missing_values_cfDNA, cbind, variants.percent = c("0"))
print(missing_values_cfDNA)

cfDNA_new <- Map(rbind, cfDNA_character, missing_values_cfDNA)
View(cfDNA_new)
print(cfDNA_new)

#combining new cfDNA and new tissue samples
SRR_117_tissue_cfDNA_new <- list(tissue_new, cfDNA_new)
names(SRR_117_tissue_cfDNA_new) <- c("tissue", "cfDNA")
print(SRR_117_tissue_cfDNA_new)

melt_117_new <- melt(SRR_117_tissue_cfDNA_new, id = "pattern", value.name = "variants.percent")
View(melt_117_new)

#converting column to a numeric one
melt_117_new$variants.percent <- as.numeric(as.character(melt_117_new$variants.percent))
str(melt_117_new)
View(melt_117_new)

#finding max in a table to put on y-axis
max_117 <- max(melt_117_new[melt_117_new$pattern %in% melt_117_filtered_variants_no_dupl$pattern,][, "variants.percent"])
print(max_117)

str(melt_117_new[melt_117_new$pattern %in% melt_117_filtered_variants_no_dupl$pattern,])
#making boxplots
testplot_2 <- ggplot(melt_117_new[melt_117_new$pattern %in% melt_117_filtered_variants_no_dupl$pattern,],aes(x = L1, y = as.numeric(variants.percent))) + 
  geom_boxplot(aes(fill=L1), outlier.shape=NA) + 
  theme(text = element_text(size = 30), legend.position = "None") +
  labs(y="Percent of variants", x="Sample origin") +
  facet_grid(. ~pattern, space = 'free_y', ) +
  coord_cartesian(ylim = c(0, as.numeric(max_117))) +
  geom_jitter() +
  ggtitle("Variants distribution for SRR_117") +
  ggeasy::easy_center_title()

ggsave('testplot_2.pdf', height = 7, width = 40)

View(melt_117_new)
p_values_117 <- melt_117_pval %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(variants.percent)) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[1]],values[[2]])$p.value)

#formatting p_values
new <- data.frame(
  group1 = "cfDNA",
  group2 = "tissue",
  y.position = 4.5
)
p_values_117 <- cbind(p_values_117, new)
print(p_values_117)

melt_117_pval <- melt_117_new[melt_117_new$pattern %in% melt_117_filtered_variants_no_dupl$pattern,]

testplot_2 + add_pvalue(p_values_117, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5)

#making boxplots for SRR_119
##variants parsing

##SRR_119
##computel_119_trimmed_no_ref
filenames_119_trimmed <- list.files("~/Desktop/tel_variants/paired/computel_119_trimmed_no_ref/output_1", pattern="SRR*", full.names=FALSE)
print(filenames_119_trimmed)
##output_1
variants_119_trimmed_1 <- list.files("~/Desktop/tel_variants/paired/computel_119_trimmed_no_ref/output_1", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
print(variants_119_trimmed_1)
trimmed_119_1 <- lapply(variants_119_trimmed_1, read.delim2)
trimmed_119_1 <- lapply(trimmed_119_1, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(trimmed_119_1)
names(trimmed_119_1) <- filenames_119_trimmed
##output_2
variants_119_trimmed_2 <- list.files("~/Desktop/tel_variants/paired/computel_119_trimmed_no_ref/output_2", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
trimmed_119_2 <- lapply(variants_119_trimmed_2, read.delim2)
trimmed_119_2 <- lapply(trimmed_119_2, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(variants_119_trimmed_2)
print(trimmed_119_2)
names(trimmed_119_2) <- filenames_119_trimmed

#combining output_1 and output_2
m <- match(names(trimmed_119_2), names(trimmed_119_1), nomatch = 0L)
print(m)
combined_119_trimmed_1_2 <- Map(rbind, trimmed_119_1[m], trimmed_119_2)
print(combined)

combined_119_trimmed_1_2 <- lapply(combined_119_trimmed_1_2, function(x) summarise_at(group_by(x, pattern), vars(abs.num), list(~floor(mean(.,na.rm=FALSE)))))
print(combined_119_trimmed_1_2)

percent_of_variants <- lapply(combined_119_trimmed_1_2, function(x) x[,sapply(x, is.numeric)]/colSums(x[,sapply(x, is.numeric)])*100)
print(percent_of_variants)
colname <- ("variants.percent")
percent_of_variants <- lapply(percent_of_variants, setNames, colname)
print(percent_of_variants)
combined_119_trimmed_1_2 <- Map(cbind, combined_119_trimmed_1_2, percent_of_variants) 
print(combined_119_trimmed_1_2)
combined_119_trimmed_1_2 <- lapply(combined_119_trimmed_1_2, function(x){
  x[order(x$abs.num, decreasing = T),]
})
print(combined_119_trimmed_1_2)
#combined_119_trimmed_1_2_top_5 <- lapply(combined_119_trimmed_1_2, function(x) x[c(2:6), ])
#combined_119_trimmed_1_2_subset <- lapply(combined_119_trimmed_1_2, function(x) subset(x, variants.percent > 1))
combined_119_trimmed_1_2_subset <- lapply(combined_119_trimmed_1_2, function(x) subset(x, pattern != "TTAGGG"))
print(combined_119_trimmed_1_2_subset)

##computel_119_untrimmed_no_ref
filenames_119_untrimmed <- list.files("~/Desktop/tel_variants/paired/computel_119_untrimmed_no_ref/output_1", pattern="SRR*", full.names=FALSE)
print(filenames_119_untrimmed)
##output_1
variants_119_untrimmed_1 <- list.files("~/Desktop/tel_variants/paired/computel_119_untrimmed_no_ref/output_1", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
print(variants_119_untrimmed_1)
untrimmed_119_1 <- lapply(variants_119_untrimmed_1, read.delim2)
untrimmed_119_1 <- lapply(untrimmed_119_1, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(untrimmed_119_1)
names(untrimmed_119_1) <- filenames_119_untrimmed
##output_2
variants_119_untrimmed_2 <- list.files("~/Desktop/tel_variants/paired/computel_119_untrimmed_no_ref/output_2", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
untrimmed_119_2 <- lapply(variants_119_untrimmed_2, read.delim2)
untrimmed_119_2 <- lapply(untrimmed_119_2, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(variants_119_untrimmed_2)
print(untrimmed_119_2)
names(untrimmed_119_2) <- filenames_119_untrimmed

#combining output_1 and output_2
q <- match(names(untrimmed_119_2), names(untrimmed_119_1), nomatch = 0L)
print(q)
combined_119_untrimmed_1_2 <- Map(rbind, untrimmed_119_1[q], untrimmed_119_2)
print(combined_119_untrimmed_1_2)

combined_119_untrimmed_1_2 <- lapply(combined_119_untrimmed_1_2, function(x) summarise_at(group_by(x, pattern), vars(abs.num), list(~floor(mean(.,na.rm=FALSE)))))
print(combined_119_untrimmed_1_2)

percent_of_variants <- lapply(combined_119_untrimmed_1_2, function(x) x[,sapply(x, is.numeric)]/colSums(x[,sapply(x, is.numeric)])*100)
print(percent_of_variants)
colname <- ("variants.percent")
percent_of_variants <- lapply(percent_of_variants, setNames, colname)
print(percent_of_variants)
combined_119_untrimmed_1_2 <- Map(cbind, combined_119_untrimmed_1_2, percent_of_variants) 
print(combined_119_untrimmed_1_2)
combined_119_untrimmed_1_2 <- lapply(combined_119_untrimmed_1_2, function(x){
  x[order(x$abs.num, decreasing = T),]
})
print(combined_119_untrimmed_1_2)
View(combined_119_untrimmed_1_2)
#combined_119_untrimmed_1_2_top_5 <- lapply(combined_119_untrimmed_1_2, function(x) x[c(2:6), ])
#combined_119_untrimmed_1_2_subset <- lapply(combined_119_untrimmed_1_2, function(x) subset(x, variants.percent > 1))
combined_119_untrimmed_1_2_subset <- lapply(combined_119_untrimmed_1_2, function(x) subset(x, pattern != "TTAGGG"))
print(combined_119_untrimmed_1_2_subset)
#combining 119 everything
variants_119 <- list(combined_119_untrimmed_1_2_subset, combined_119_trimmed_1_2_subset)
variants_119 <- do.call(c, variants_119)
variants_119 <- lapply(variants_119, function(x) x[(names(x) %in% c("pattern", "variants.percent"))])
View(variants_119)
print(variants_119)

##combining dataframes from cancer/healthy, tissue/cfDNA within datasets

##SRR_119 -- cancer/healthy
cancer <- list_merge(variants_119[c(2:6)]) 
View(cancer)
healthy <- list_merge(variants_119[1])
View(healthy)
SRR_119_cancer_healthy <- list(cancer, healthy)
names(SRR_119_cancer_healthy) <- c("cancer", "healthy")
View(SRR_119_cancer_healthy)

##melting the data and making plots
##SRR_119
melt_119 <- melt(SRR_119_cancer_healthy, id = "pattern", value.name = "variants.percent") 
View(melt_119)

#should I filter before or after grouping by?
#filtering before grouping by

melt_119_filtered <- subset(melt_119, variants.percent > 1)
View(melt_119_filtered)

#grouping by repeated variants to get a table with possible variants (>1%)
melt_119_filtered_variants <- melt_119_filtered %>% group_by(pattern, L1) %>%
  summarise(variants.percent=(mean(as.numeric(variants.percent))))
View(melt_119_filtered_variants)

#dealing with duplicates, this way you would just have variants with a percent in samples > 1%
melt_119_filtered_variants_no_dupl <- distinct(melt_119_filtered_variants, pattern)
View(melt_119_filtered_variants_no_dupl)

#calculating the average over all the samples, including ones with a percent < 1
melt_119_average <- melt_119[melt_119$pattern %in% melt_119_filtered_variants_no_dupl$pattern,] %>% group_by(pattern, L1) %>% summarise(variants.percent=(mean(as.numeric(variants.percent))))
View(melt_119_average)

#a different approach with zeros
###the most common variants for cancer/healthy, tissue/cfDNA

#SRR_119

#there is a need to put zeros in variants percent if the certain variant doesn't exist in a sample
#if there is no variant from the melt_117_grouped for a sample, then set it to zero
#cancer samples
missing_values_cancer <- lapply(cancer, function(x) anti_join(melt_119_filtered_variants_no_dupl[,"pattern"], x, by = 'pattern'))
missing_values_cancer <- lapply(missing_values_cancer, cbind, variants.percent = c("0"))
print(missing_values_cancer)

cancer_new <- Map(rbind, cancer, missing_values_cancer)
View(cancer_new)
print(cancer_new)

#healthy samples
healthy_character <- outer(healthy)
print(healthy_character)

missing_values_healthy <- lapply(healthy_character, function(x) anti_join(melt_119_filtered_variants_no_dupl[,"pattern"], x, by = 'pattern'))

missing_values_healthy <- lapply(missing_values_healthy, cbind, variants.percent = c("0"))
print(missing_values_healthy)

healthy_new <- Map(rbind, healthy_character, missing_values_healthy)
View(healthy_new)
print(healthy_new)

#combining new cfDNA and new tissue samples
SRR_119_cancer_healthy_new <- list(cancer_new, healthy_new)
names(SRR_119_cancer_healthy_new) <- c("cancer", "healthy")
print(SRR_119_cancer_healthy_new)
str(SRR_119_cancer_healthy_new)

melt_119_new <- melt(SRR_119_cancer_healthy_new, id = "pattern", value.name = "variants.percent")
View(melt_119_new)

#converting column to a numeric
melt_119_new$variants.percent <- as.numeric(as.character(melt_119_new$variants.percent))
str(melt_119_new)
View(melt_119_new)

#finding max in a table to put on y-axis
max_119 <- max(melt_119_new[melt_119_new$pattern %in% melt_119_filtered_variants_no_dupl$pattern,][, "variants.percent"])
print(max_119)

#making boxplots
testplot_3 <- ggplot(melt_119_new[melt_119_new$pattern %in% melt_119_filtered_variants_no_dupl$pattern,],aes(x = L1, y = as.numeric(variants.percent))) + 
  geom_boxplot(aes(fill=L1), outlier.shape=NA) + 
  theme(text = element_text(size = 30), legend.position = "None") +
  labs(y="Percent of variants", x="Disease status") +
  facet_grid(. ~pattern, space = 'free_y', ) +
  coord_cartesian(ylim = c(0, as.numeric(max_119))) +
  geom_jitter() +
  ggtitle("Variants distribution for SRR_119") +
  ggeasy::easy_center_title()

ggsave('testplot_3.pdf', height = 7, width = 60, limitsize = FALSE)

melt_119_pval <- melt_119_new[melt_119_new$pattern %in% melt_119_filtered_variants_no_dupl$pattern,]
print(melt_119_pval)

p_values_119 <- melt_119_pval %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(variants.percent)) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[1]],values[[2]])$p.value)

#formatting p_values
new <- data.frame(
  group1 = "cancer",
  group2 = "healthy",
  y.position = 3.5
)
p_values_119 <- cbind(p_values_119, new)
print(p_values_119)

testplot_3 + add_pvalue(p_values_119, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5)
ggsave('testplot_3.pdf', height = 7, width = 60, limitsize = FALSE)

#making boxplots for SRR_744
##SRR_744
##computel_744_trimmed_no_ref
filenames_744_trimmed <- list.files("~/Desktop/tel_variants/paired/computel_744_trimmed_no_ref/output_1", pattern="SRR*", full.names=FALSE)
print(filenames_744_trimmed)
##output_1
variants_744_trimmed_1 <- list.files("~/Desktop/tel_variants/paired/computel_744_trimmed_no_ref/output_1", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
print(variants_744_trimmed_1)
trimmed_744_1 <- lapply(variants_744_trimmed_1, read.delim2)
trimmed_744_1 <- lapply(trimmed_744_1, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(trimmed_744_1)
names(trimmed_744_1) <- filenames_744_trimmed
##output_2
variants_744_trimmed_2 <- list.files("~/Desktop/tel_variants/paired/computel_744_trimmed_no_ref/output_2", pattern="tel.variants.txt", full.names=TRUE, recursive = TRUE)
trimmed_744_2 <- lapply(variants_744_trimmed_2, read.delim2)
trimmed_744_2 <- lapply(trimmed_744_2, function(x) x[(names(x) %in% c("pattern", "abs.num", "X..of.all.patterns"))])
print(variants_744_trimmed_2)
print(trimmed_744_2)
names(trimmed_744_2) <- filenames_744_trimmed

#combining output_1 and output_2
m <- match(names(trimmed_744_2), names(trimmed_744_1), nomatch = 0L)
print(m)
combined_744_trimmed_1_2 <- Map(rbind, trimmed_744_1[m], trimmed_744_2)
print(combined)

combined_744_trimmed_1_2 <- lapply(combined_744_trimmed_1_2, function(x) summarise_at(group_by(x, pattern), vars(abs.num), list(~floor(mean(.,na.rm=FALSE)))))
df <- lapply(combined_744_trimmed_1_2, function(x) as.data.frame(x))
print(df)
print(combined_744_trimmed_1_2)

percent_of_variants <- lapply(combined_744_trimmed_1_2, function(x) x[,sapply(x, is.numeric)]/colSums(x[,sapply(x, is.numeric)])*100)
print(percent_of_variants)
colname <- ("variants.percent")
percent_of_variants <- lapply(percent_of_variants, setNames, colname)
print(percent_of_variants)
combined_744_trimmed_1_2 <- Map(cbind, combined_744_trimmed_1_2, percent_of_variants) 
print(combined_744_trimmed_1_2)
combined_744_trimmed_1_2 <- lapply(combined_744_trimmed_1_2, function(x){
  x[order(x$abs.num, decreasing = T),]
})
print(combined_744_trimmed_1_2)
#variants_744 <- lapply(combined_744_trimmed_1_2, function(x) x[c(2:6), ])
#variants_744 <- lapply(combined_744_trimmed_1_2, function(x) subset(x, variants.percent > 1))
variants_744 <- lapply(combined_744_trimmed_1_2, function(x) subset(x, pattern != "TTAGGG"))
print(variants_744)
View(variants_744)
variants_744 <- lapply(variants_744, function(x) x[(names(x) %in% c("pattern", "variants.percent"))])

##SRR_744 -- cancer pre-, post-operation/healthy
healthy <- list_merge(variants_744[c(17:20)])
View(healthy)
cancer_preop <- list_merge(variants_744[c(1:11)])
cancer_postop <- list_merge(variants_744[c(12:16)])
View(cancer_preop)
View(cancer_postop)
SRR_744_cancer_healthy <- list(cancer_preop, cancer_postop, healthy)
View(SRR_744_cancer_healthy)
names(SRR_744_cancer_healthy) <- c("cancer_preop", "cancer_postop", "healthy")

melt_744 <- melt(SRR_744_cancer_healthy, id = "pattern", value.name = "variants.percent")
View(melt_744)

#should I filter before or after grouping by?
#filtering before grouping by

melt_744_filtered <- subset(melt_744, variants.percent > 1)
View(melt_744_filtered)

#grouping by repeated variants to get a table with possible variants (>1%)
melt_744_filtered_variants <- melt_744_filtered %>% group_by(pattern, L1) %>%
  summarise(variants.percent=(mean(as.numeric(variants.percent))))
View(melt_744_filtered_variants)

#dealing with duplicates, this way you would just have variants with a percent in samples > 1%
melt_744_filtered_variants_no_dupl <- distinct(melt_744_filtered_variants, pattern)
View(melt_744_filtered_variants_no_dupl)

#calculating the average over all the samples, including ones with a percent < 1
melt_744_average <- melt_744[melt_744$pattern %in% melt_744_filtered_variants_no_dupl$pattern,] %>% group_by(pattern, L1) %>% summarise(variants.percent=(mean(as.numeric(variants.percent))))
View(melt_744_average)

#a different approach with zeros
###the most common variants for cancer/healthy, tissue/cfDNA

#SRR_744

#there is a need to put zeros in variants percent if the certain variant doesn't exist in a sample
#if there is no variant from the melt_117_grouped for a sample, then set it to zero
#cancer_preop samples
missing_values_cancer_preop <- lapply(cancer_preop, function(x) anti_join(melt_744_filtered_variants_no_dupl[,"pattern"], x, by = 'pattern'))
missing_values_cancer_preop <- lapply(missing_values_cancer_preop, cbind, variants.percent = c("0"))
print(missing_values_cancer_preop)

cancer_preop_new <- Map(rbind, cancer_preop, missing_values_cancer_preop)
View(cancer_preop_new)
print(cancer_preop_new)

missing_values_cancer_postop <- lapply(cancer_postop, function(x) anti_join(melt_744_filtered_variants_no_dupl[,"pattern"], x, by = 'pattern'))

missing_values_cancer_postop <- lapply(missing_values_cancer_postop, cbind, variants.percent = c("0"))
print(missing_values_healthy)

cancer_postop_new <- Map(rbind, cancer_postop, missing_values_cancer_postop)
View(cancer_postop_new)
print(cancer_postop_new)

#healthy
missing_values_healthy <- lapply(healthy, function(x) anti_join(melt_744_filtered_variants_no_dupl[,"pattern"], x, by = 'pattern'))
missing_values_healthy <- lapply(missing_values_healthy, cbind, variants.percent = c("0"))
print(missing_values_healthy)

healthy_new <- Map(rbind, healthy, missing_values_healthy)
View(healthy_new)
print(healthy_new)

#combining new cfDNA and new tissue samples
SRR_744_cancer_healthy_new <- list(cancer_preop_new, cancer_postop_new, healthy_new)
names(SRR_744_cancer_healthy_new) <- c("cancer_preop", "cancer_postop", "healthy")
print(SRR_744_cancer_healthy_new)
str(SRR_744_cancer_healthy_new)

melt_744_new <- melt(SRR_744_cancer_healthy_new, id = "pattern", value.name = "variants.percent")
View(melt_744_new)

#converting column to a numeric
melt_744_new$variants.percent <- as.numeric(as.character(melt_744_new$variants.percent))
str(melt_744_new)
View(melt_744_new)

#finding max in a table to put on y-axis
max_744 <- max(melt_744_new[melt_744_new$pattern %in% melt_744_filtered_variants_no_dupl$pattern,][,"variants.percent"])
print(max_744)

#making boxplots
testplot_4 <- ggplot(melt_744_new[melt_744_new$pattern %in% melt_744_filtered_variants_no_dupl$pattern,],aes(x = L1, y = as.numeric(variants.percent))) +
  geom_boxplot(aes(fill=L1), outlier.shape=NA) +
  theme(text = element_text(size = 30), legend.position = "None") +
  labs(y="Percent of variants", x="Disease status") +
  facet_grid(. ~pattern, space = 'free_y', ) +
  coord_cartesian(ylim = c(0, as.numeric(max_744))) +
  geom_jitter() +
  ggtitle("Variants distribution for SRR_744") +
  ggeasy::easy_center_title()

ggsave('testplot_4.pdf', height = 10, width = 100, limitsize = FALSE)

melt_744_pval <- melt_744_new[melt_744_new$pattern %in% melt_744_filtered_variants_no_dupl$pattern,]
print(melt_744_pval)
View(melt_744_pval)

str(melt_744_pval)
p_values_744_1 <- melt_744_pval %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(variants.percent)) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[1]],values[[2]])$p.value)

p_values_744_2 <- melt_744_pval %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(variants.percent)) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[2]],values[[3]])$p.value)

p_values_744_3 <- melt_744_pval %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(variants.percent)) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[1]],values[[3]])$p.value)

#formatting p_values
new_1 <- data.frame(
  group1 = "cancer_postop",
  group2 = "cancer_preop",
  y.position = 10
)
p_values_744_1 <- cbind(p_values_744_1, new_1)
print(p_values_744_1)

new_2 <- data.frame(
  group1 = "cancer_preop",
  group2 = "healthy",
  y.position = 10.5
)

p_values_744_2 <- cbind(p_values_744_2, new_2)
print(p_values_744_2)

new_3 <- data.frame(
  group1 = "cancer_postop",
  group2 = "healthy",
  y.position = 11.5
)

p_values_744_3 <- cbind(p_values_744_3, new_3)
print(p_values_744_3)

testplot_4 + add_pvalue(p_values_744_1, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5) + 
  add_pvalue(p_values_744_2, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5) +
  add_pvalue(p_values_744_3, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5)

ggsave('testplot_4.pdf', height = 10, width = 100, limitsize = FALSE)

#making boxplots for canonical variants
#SRR_117
single_subset_can <- lapply(single, function(x) subset(x, pattern == "TTAGGG"))
column_names <- c("pattern", "abs.num", "variants.percent")
single_subset_can <- lapply(single_subset_can, setNames, column_names)
names(single_subset_can) <- filenames_single
print(single_subset_can)

combined_117_trimmed_1_2_subset_can <- lapply(combined_117_trimmed_1_2, function(x) subset(x, pattern == "TTAGGG"))
print(combined_117_trimmed_1_2_subset_can)

combined_117_untrimmed_1_2_subset_can <- lapply(combined_117_untrimmed_1_2, function(x) subset(x, pattern == "TTAGGG"))
print(combined_117_untrimmed_1_2_subset_can)

variants_117_can <- list(combined_117_untrimmed_1_2_subset_can, combined_117_trimmed_1_2_subset_can, single_subset_can)
View(variants_117)
print(variants_117_can)
variants_117_can <- do.call(c, variants_117_can)
variants_117_can <- lapply(variants_117_can, function(x) x[(names(x) %in% c("pattern", "variants.percent"))])
View(variants_117_can)

tissue <- list_merge(variants_117_can[c(1, 19, 5, 20, 21, 22, 23, 24)])
cfDNA <- list_merge(variants_117_can[-c(1, 19, 5, 20, 21, 22, 23, 24)])
View(cfDNA)
SRR_117_cancer_healthy <- list(tissue, cfDNA)
names(SRR_117_cancer_healthy) <- c("tissue", "cfDNA")
print(SRR_117_cancer_healthy)

melt_117 <- melt(SRR_117_cancer_healthy, id = "pattern", value.name = "variants.percent")
View(melt_117)

#finding max in a table to put on y-axis
max_117 <- max(melt_117[, "variants.percent"])
min_117 <- min(melt_117[, "variants.percent"])
print(min_117)

#making boxplots
testplot_5 <- ggplot(melt_117, aes(x = L1, y = as.numeric(variants.percent))) +
  geom_boxplot(aes(fill=L1), outlier.shape=NA) +
  theme(text = element_text(size = 30), legend.position = "None") +
  labs(y="Percent of variants", x="Sample origin") +
  coord_cartesian(ylim = c(as.numeric(min_117), as.numeric(max_117))) +
  geom_jitter() +
  ggtitle("TTAGGG distribution for SRR_117") +
  ggeasy::easy_center_title()

ggsave('testplot_5.pdf', height = 5, width = 10)
str(melt_117)
p_values_117 <- melt_117 %>%
  group_by(pattern,L1) %>%
  summarise(values = list(as.numeric(variants.percent))) %>%
  group_by(pattern) %>%
  summarise(p_value= wilcox.test(values[[1]],values[[2]])$p.value)
print(p_values_117)
#formatting p_values
new <- data.frame(
  group1 = "cfDNA",
  group2 = "tissue",
  y.position = 95
)
p_values_117 <- cbind(p_values_117, new)
print(p_values_117)

testplot_5 + add_pvalue(p_values_117, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5)

#SRR_119
combined_119_trimmed_1_2_subset <- lapply(combined_119_trimmed_1_2, function(x) subset(x, pattern == "TTAGGG"))
combined_119_untrimmed_1_2_subset <- lapply(combined_119_untrimmed_1_2, function(x) subset(x, pattern == "TTAGGG"))
print(combined_119_trimmed_1_2_subset)
variants_119 <- list(combined_119_untrimmed_1_2_subset, combined_119_trimmed_1_2_subset)
variants_119 <- do.call(c, variants_119)
variants_119 <- lapply(variants_119, function(x) x[(names(x) %in% c("pattern", "variants.percent"))])
View(variants_119)
print(variants_119)

##combining dataframes from cancer/healthy, tissue/cfDNA within datasets

##SRR_119 -- cancer/healthy
cancer <- list_merge(variants_119[c(2:6)]) 
View(cancer)
healthy <- list_merge(variants_119[1])
View(healthy)
SRR_119_cancer_healthy <- list(cancer, healthy)
names(SRR_119_cancer_healthy) <- c("cancer", "healthy")
View(SRR_119_cancer_healthy)

##melting the data and making plots
##SRR_119
melt_119 <- melt(SRR_119_cancer_healthy, id = "pattern", value.name = "variants.percent") 
View(melt_119)

#finding max in a table to put on y-axis
max_119 <- max(melt_119[, "variants.percent"])
min_119 <- min(melt_119[, "variants.percent"])
print(min_119)
print(max_119)

#making boxplots
testplot_6 <- ggplot(melt_119,aes(x = L1, y = as.numeric(variants.percent))) + 
  geom_boxplot(aes(fill=L1), outlier.shape=NA) + 
  theme(text = element_text(size = 30), legend.position = "None") +
  labs(y="Percent of variants", x="Disease status") +
  coord_cartesian(ylim = c(as.numeric(min_119), as.numeric(max_119))) +
  geom_jitter() +
  ggtitle("TTAGGG distribution for SRR_119") +
  ggeasy::easy_center_title()

ggsave('testplot_6.pdf', height = 5, width = 10)

p_values_119 <- melt_119 %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(as.numeric(variants.percent))) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[1]],values[[2]])$p.value)
print(p_values_119)
#formatting p_values
new <- data.frame(
  group1 = "cancer",
  group2 = "healthy",
  y.position = 85
)
p_values_119 <- cbind(p_values_119, new)
print(p_values_119)

testplot_6 + add_pvalue(p_values_119, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5)
ggsave('testplot_3.pdf', height = 7, width = 60, limitsize = FALSE)

#SRR_744
variants_744 <- lapply(combined_744_trimmed_1_2, function(x) subset(x, pattern == "TTAGGG"))
print(variants_744)
View(variants_744)
variants_744 <- lapply(variants_744, function(x) x[(names(x) %in% c("pattern", "variants.percent"))])

##SRR_744 -- cancer pre-, post-operation/healthy
healthy <- list_merge(variants_744[c(17:20)])
View(healthy)
cancer_preop <- list_merge(variants_744[c(1:11)])
cancer_postop <- list_merge(variants_744[c(12:16)])
View(cancer_preop)
View(cancer_postop)
SRR_744_cancer_healthy <- list(cancer_preop, cancer_postop, healthy)
View(SRR_744_cancer_healthy)
names(SRR_744_cancer_healthy) <- c("cancer_preop", "cancer_postop", "healthy")

melt_744 <- melt(SRR_744_cancer_healthy, id = "pattern", value.name = "variants.percent")
View(melt_744)

#finding max in a table to put on y-axis
max_744 <- max(melt_744[, "variants.percent"])
min_744 <- min(melt_744[, "variants.percent"])
print(min_744)
print(max_744)

#making boxplots
testplot_7 <- ggplot(melt_744,aes(x = L1, y = as.numeric(variants.percent))) + 
  geom_boxplot(aes(fill=L1), outlier.shape=NA) + 
  theme(text = element_text(size = 30), legend.position = "None") +
  labs(y="Percent of variants", x="Disease status") +
  coord_cartesian(ylim = c(as.numeric(min_744), as.numeric(max_744))) +
  geom_jitter() +
  ggtitle("TTAGGG distribution for SRR_744") +
  ggeasy::easy_center_title()

ggsave('testplot_7.pdf', height = 5, width = 10)

p_values_744_1 <- melt_744 %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(as.numeric(variants.percent))) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[1]],values[[2]])$p.value)

p_values_744_2 <- melt_744 %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(as.numeric(variants.percent))) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[2]],values[[3]])$p.value)

p_values_744_3 <- melt_744 %>% 
  group_by(pattern,L1) %>% 
  summarise(values = list(as.numeric(variants.percent))) %>%
  group_by(pattern) %>% 
  summarise(p_value= wilcox.test(values[[1]],values[[3]])$p.value)

#formatting p_values
new_1 <- data.frame(
  group1 = "cancer_postop",
  group2 = "cancer_preop",
  y.position = 90
)
p_values_744_1 <- cbind(p_values_744_1, new_1)
print(p_values_744_1)

new_2 <- data.frame(
  group1 = "cancer_preop",
  group2 = "healthy",
  y.position = 93
)

p_values_744_2 <- cbind(p_values_744_2, new_2)
print(p_values_744_2)

new_3 <- data.frame(
  group1 = "cancer_postop",
  group2 = "healthy",
  y.position = 95
)

p_values_744_3 <- cbind(p_values_744_3, new_3)
print(p_values_744_3)

testplot_7 + add_pvalue(p_values_744_1, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5) + 
  add_pvalue(p_values_744_2, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5) +
  add_pvalue(p_values_744_3, label = "p = {round(p_value, digits=5)}", remove.bracket = TRUE, label.size =6.5)

ggsave('testplot_7.pdf', height = 5, width = 10)
