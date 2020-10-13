library(here)
oldName = list.files(here("data/trimmed"),
                     pattern = "R1.fastq.trim.gz",
                     full.names = TRUE)
newName = gsub("_R1.fastq.trim.gz",".trim.R1.fastq.gz", oldName)
file.rename(oldName, newName)

oldName = list.files(here("data/trimmed"),
                     pattern = "R2.fastq.trim.gz",
                     full.names = TRUE)
newName = gsub("_R2.fastq.trim.gz",".trim.R2.fastq.gz", oldName)
file.rename(oldName, newName)

oldName = list.files(here("data/trimmed"),
                     pattern = "R1.fastq.trim.unpaired.gz",
                     full.names = TRUE)
newName = gsub("_R1.fastq.trim.unpaired.gz",".trim.unpaired.R1.fastq.gz", oldName)
file.rename(oldName, newName)

oldName = list.files(here("data/trimmed"),
                     pattern = "R2.fastq.trim.unpaired.gz",
                     full.names = TRUE)
newName = gsub("_R2.fastq.trim.unpaired.gz",".trim.unpaired.R2.fastq.gz", oldName)
file.rename(oldName, newName)
