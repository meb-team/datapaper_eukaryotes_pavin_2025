library(data.table)
library(tidyverse)

# Aim of the script :
# - Get a number of mapped read per MAG
# - Compute the detection (breadth of coverage) for a MAG

# 0. Setup
set.seed(6843651)

# 1.1 Read the Contig2MAG table ----
## A two-columns tab-delimited file with the contig name and its MAG
## Can be very long if you have a lots of MAG
ctg2mag <- fread("contigs_to_MAGs.tsv",
  sep = "\t", header = FALSE,
  col.names = c("contig", "MAG")
)

# # 1.2 Read the Size2MAG table ----
mag_length <- fread("MAGs_length.tsv",
  sep = "\t", header = TRUE
)

# 2. Read the result ----
## UPDATE THE "pattern" ARGUMENT TO MATCH ALL YOUR FILES
indata <- "mapping/" # Directory where to search all files
files <- list.files(indata,
  pattern = "*.filter.cov.tsv.gz",
  full.names = TRUE, recursive = TRUE
)

# 3. Main loop ----
## Keep the number of read mapped per MAG/SAG
## I would also use the detection (SAMtools 'coverage' column)

first <- TRUE
count <- 0

### Loop to get the detection per MAG per Metagenome ----
for (file in files) {
  sample <- unlist(strsplit(basename(file), "[.]"))[1]

  # read
  dt <- fread(file)

  # Check files is not empty, otherwise Next iteration
  if (dim(dt)[1] == 0) {
    print(paste0("    empty sample: ", sample))
    next
  }

  # Add the MAG info, aka merge on "contig"
  dt <- dt[ctg2mag, on = .("#rname" = contig)]

  # Compute the detection:
  ## sum of covered bases / MAG length then multiply by 100
  dt_detec <- dt[, .(detection = (sum(covbases) / sum(endpos)) * 100), by = .(MAG)]
  setnames(dt_detec, "detection", sample)

  # Merge the result
  if (first == TRUE) {
    first <- FALSE
    res_detec <- copy(dt_detec)
  } else {
    res_detec <- res_detec[dt_detec, on = .(MAG = MAG)]
  }

  # Clean
  rm(list = c("dt", "dt_detec"))

  # Log
  count <- count + 1
  if (count %% 100 == 0) {
    print(paste0("Processed ", as.character(count), " files"))
  }
}
# Save the table
fwrite(res_detec,
  "detection_per_MAG.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)

## Loop that only sum the number of reads ----

for (file in files) {
  sample <- unlist(strsplit(basename(file), "[.]"))[1]

  # read
  dt <- fread(file)

  # Check files is not empty, otherwise Next iteration
  if (dim(dt)[1] == 0) {
    print(paste0("    empty sample: ", sample))
    next
  }

  # Add the MAG info, aka merge on "contig"
  dt_mag <- dt[ctg2mag, on = .("#rname" = contig)]

  # Compute the sum of read per MAG
  dt_mag_sum <- dt_mag[, sum(numreads), by = MAG]
  setnames(dt_mag_sum, "V1", sample)

  # Merge the result
  if (first == TRUE) {
    first <- FALSE
    res <- copy(dt_mag_sum)
  } else {
    res <- res[dt_mag_sum, on = .(MAG = MAG)]
  }

  # Clean
  rm(list = c("dt_mag", "dt_mag_sum"))

  # Log
  count <- count + 1
  if (count %% 100 == 0) {
    print(paste0("Processed ", as.character(count), " files"))
  }
}

# Save the table
fwrite(res,
  "sum_read_per_MAG.tsv",
  sep = "\t", quote = FALSE, row.names = FALSE
)
