# Data-preparation-wrangling
For scripts to run data preparation, wrangling and processing. All codes written in R (so far).

## Data wrangling using v2 pipeline data

When making a meal in the kitchen, you need 3 things: cutlery, ingredients, and a single kitchen. 

Likewise, before even cooking (or running your analysis), you need your cutlery (R packages), ingredients (your data), and your cooking station (working directory).

First, let's (prepare the cutlery) load the required packages for microbiome analysis.

```{r cars}
library(tidyverse)
library(vegan)
library(phyloseq)
library(Maaslin2)
library(ggsci)
library(microbiome)
library(data.table)
library(ggrepel)

# Then, we (pick a kitchen to cook in, where your ingredients are stored) set your working directory, where your data is stored.
setwd("/Users/IanLim/foldername")
current_dir <- getwd() # Get the current working directory
list.files() # Let's make sure our data is there
```

For metagenomic analyses, we need to prepare a single object that allows us to seamlessly run whatever analysis step we need, without linking multiple files such as demographic data or count data files. Everything will be in one single object, called a 'phyloseq'.

This phyloseq object consists of 3 files, which is typical in every metagenomic analysis.
1. asv: .tsv file(s) consisting of the read count data/abundance of each bacterial taxa with its sequences (rows), found in each sample (columns).
2. tax: Derived from the asv file, it consists of taxonomic information (columns: Kingdom, Phylum, Class, Order, Family, Genus, Species) of the bacterial taxa with its sequences (rows). 
3. met: Typically a .csv file, consisting of the sample IDs, demographic information, etc. of your study.

So, now that we have the 3 ingredients we need to make a phyloseq, let's prepare them.
Let's start with our first object, called 'asv'. Depending on your data, you may receive the data in batches, or in a single batch after sequencing. Today we are going to prepare it as if we have received our data in five batches.

```{r cars}
### ASV ###

# Select asv (.tsv) files in your directory/folder
asvfilestoprocess = dir(pattern = "\\.tsv$")
asvfilestoprocess # let's see if they have been correctly selected

# Import asv (.tsv) files into your directory/folder
asvfiles = lapply(asvfilestoprocess, function(x) {
  read.delim(x, header = T)
})

# Rename the components in asvfiles based on their respective asv batches
batchname = c("1", "2", "3", "4", "5")
names(asvfiles) = batchname

# Let's one of the batches 
View(asvfiles[["1"]])

# Select only asv abundance columns from asvfiles and assign to asv
asv = lapply(asvfiles, function(x) {
  x %>% select(-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")) # removing taxonomic columns, because we need them in TAX only
})

# Now, to make a single ASV file, merge all components in asv variable and assign first column (the ASV sequences) as rownames
asv = plyr::join_all(asv, by="X.seq", type="full") %>% mutate(across(.fns = ~ifelse(is.na(.), 0, .))) %>% column_to_rownames("X.seq")

# A bit of tidying up, remove “_F_filt.fastq.gz" from the Analysis_ID names in the asv columns
colnames(asv) = gsub("\\_F_filt.fastq.gz", "", colnames(asv))

```

Now that we have an ASV file, 2 more to go. 

We can use the asvfiles earlier to create a taxonomic object, let's call it 'tax'.

```{cars}
### TAX ###

# Select the taxonomic ranks (should be the last 7 columns) from asvfiles and assign them to 'tax'
tax = lapply(asvfiles, function(x) {
  x[,c(1, (ncol(x)-6):ncol(x))]
})

# Merge all components in tax variable (merging all tax tables from the many batches)
tax = plyr::join_all(tax, by = "X.seq", type= "full") %>% column_to_rownames("X.seq")

```

We have created 2/3 of our phyloseq objects.

Now, let's import our metadata file from our folder.
```{cars}
#### MET ####
met <- read.csv('metadata.csv')

# Assuming your sample IDs are under a column named AnalysisID, set AnalysisID as the same format as the AnalysisID in asv column names
# For example, in the metadata file, each AnalysisID has '-' in between each 4 letters/digits.
met$Analysis_ID <- gsub("-", "", met$Analysis_ID)

# To match the asv and met files together, set the AnalysisID column as the rownames
rownames(met) <- met$Analysis_ID
met$Analysis_ID <- NULL  # Remove the AnalysisID column to avoid redundancy
```

We are close to create our phyloseq object & begin our analysis, but first.

We need to make sure all 3 files (asv, tax, met) match each other. 

```{cars}
## Verify rownames are matching for all three files
# Check asv and tax rownames
all(rownames(asv) %in% rownames(tax))

# Check asv columnnames and met rownames
all(colnames(asv) %in% rownames(met)) # It will appear false, because not all the samples in the batches are needed

# Filter the ASV data to include only the sample IDs in the list
sample_ids_to_keep <- rownames(met) # List of 16 sample IDs from the "met" dataframe
filtered_asv <- asv[, sample_ids_to_keep] # Filter

# Now, the "filtered_asv" dataframe contains only the 16 sample IDs you need.
all(colnames(filtered_asv) %in% rownames(met))
```

When we analyze our bacterial data, and identify specific ASVs, there would be multiple ASVs of the same species. Hence, let's make sure each of them are individually identifiable. We will name each ASV row with the highest level of taxonomic annotation detected (eg. Species level is the lowest & most specific annotation possible).

```{cars}
# Rename row names and attach respective taxa names
attach_taxa_name <- function(row, taxa) {
  if (taxa == "Species") {
    if (!is.na(tax[row, "Species"])) {
      return(paste0("ASV", row, "_", tax[row, "Genus"], "_", tax[row, "Species"]))
    } else {
      return(paste0("ASV", row, "_", tax[row, "Genus"]))
    }
  } else {
    return(paste0("ASV", row, "_", tax[row, taxa]))
  }
}

# Get the lowest level of annotation for each row
highest_taxa <- apply(tax[, c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")], 1, function(row) {
  non_na_indices <- which(!is.na(row))
  if (length(non_na_indices) == 0) {
    return("Unknown")
  } else {
    max_idx <- max(non_na_indices)
    if (max_idx == 7) {
      return("Species")
    } else {
      return(names(row)[max_idx])
    }
  }
})

# Rename row names and attach the respective taxa names
row_names <- mapply(attach_taxa_name, seq_len(nrow(tax)), highest_taxa)
rownames(filtered_asv) <- ifelse(is.na(row_names), rownames(filtered_asv), row_names)
rownames(tax) <- ifelse(is.na(row_names), rownames(tax), row_names)

```

Let's create our phyloseq object!

```{cars}
## Now with the 3 main ingredients, we can create our phyloseq object
asv.ps = otu_table(filtered_asv, taxa_are_rows=T)
tax.ps = tax_table(as.matrix(tax))
met.ps = sample_data(met)

ps <- phyloseq(asv.ps, tax.ps, met.ps)
```
