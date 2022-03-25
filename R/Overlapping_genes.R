library(tidyverse)

similarity_blast <- read.csv("~/Downloads/blast_results_1-ResultsTable-Homo_sapiens_Tools_Blast_.csv")

similarity_blast$Length <- as.integer(gsub(" \\[Sequence\\]","",similarity_blast$Length))
similarity_blast$Genomic.Location <- gsub(" \\[Sequence\\]","",similarity_blast$Genomic.Location)
similarity_blast$X.ID <- as.double(gsub(" \\[Alignment\\]","",similarity_blast$X.ID))

similarity_blast <- similarity_blast %>%
  filter(!(Overlapping.Gene.s. == "")) %>%
  separate_rows(Overlapping.Gene.s.,
                sep = ", ") %>%
  separate(Genomic.Location,
           sep = ":",
           into = c("chr", "pos"))

overlapping_genes <- c(similarity_blast$Overlapping.Gene.s.) %>%
  unique()
overlapping_genes <- as_tibble(overlapping_genes[!(overlapping_genes == "APOO")])
write.table(overlapping_genes,
            file = 'overlapping_genes.txt',
            quote = FALSE, 
            sep= '\t' ,
            col.names = FALSE,
            row.names = FALSE)

converted_genes <- read.csv("~/Desktop/GSE124439_RAW/gProfiler_hsapiens_14-03-2022_14-58-51.csv")
converted_genes <- converted_genes[-72,]
write.table(converted_genes$converted_alias,
            file = 'converted_genes.txt',
            quote = FALSE, 
            sep= '\t' ,
            col.names = FALSE,
            row.names = FALSE)
