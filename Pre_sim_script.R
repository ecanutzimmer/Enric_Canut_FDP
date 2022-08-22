# Script to prepare the simulation 

library(polyester)
library(Biostrings)
library(dplyr)
library(tidyverse)
library(stringr)

# load into an object the featurecounts ouput, get rid of the unwanted columns 

counts_df <- read.table(file = "data/featureCounts.txt", header = T) %>% dplyr::select(-Chr, -Start, -End, -Strand, -Length)

# rename counts column and check dataframe

names(counts_df)[2] <- 'counts'
nrow(counts_df) 
head(counts_df)

# assign a unique identifier to each copy of a repeat. I will achieve so by adding a "_dupX" at the end of each repeat name, where X corresponds to the nth copy of an specific repeat

counts_df_2 <- counts_df %>% group_by(Geneid) %>% mutate(Geneid=paste(Geneid, row_number(), sep = "_dup")) %>% ungroup() %>% column_to_rownames("Geneid")

# sort dataframe by aplhabetical order 

rownames_df <-  cbind(row.names(counts_df_2),counts_df_2)
sorted_rownames <- rownames_df[order(rownames_df[,1]),]
names(sorted_rownames)[1] <- 'rownames' 
sorted_counts <- sorted_rownames %>% dplyr::select(-rownames)

# we convert counts_df_2 to an integer vector (counts_vect) and we sum +1 to each reads_per_transcript value since polyester has a limitation when a transcript has 0 reads associated. 

counts_vect <- as.vector(as.integer(sorted_counts[,1]+1))
length(counts_vect)
is.null(counts_vect)

# we selected 20 TE's to be known as active in the human genome and we assigned a fold change to represent over expression in one condition or the other.  

# command to check the counts associated to each copy of a repeat family 

counts_df_2 %>% tibble::rownames_to_column("Geneid") %>% dplyr::mutate(row = row_number()) %>% dplyr::filter(stringr::str_detect(Geneid, "MER20"))

# defining fold_change matrix

fc <- matrix(rep(1, 2*335744), nrow = 335744)

# first 10 transcripts to be overexpressed in group 1

fc[232554,1]=1.5  # 397 counts | MER20DNA/hAT-Charlie
fc[101887,1]=4 # 4 counts | Charlie18aDNA/hAT-Charlie
fc[141240,1]=2  # 3173 counts | L1MB2LINE/L1
fc[26663,1]=4 # 796 counts | AluJbSINE/Alu
fc[281576,1]=2  # 239 counts | MIRbSINE/MIR
fc[313119,1]=1.5  #11 counts | MLT1ILTR/ERVL-MaLR
fc[249649,1]=2  #23 counts | MER8DNA/TcMar-Tigger
fc[49727,1]=1.5  # 75 counts | AluSx1SINE/Alu
fc[119113,1]=3  # 162 counts | HERV4_ILTR/ERV1
fc[309333,1]=4 # 60 counts  | MLT1DLTR/ERVL-MaLR

# last 10 transcripts to be overxpressed in group 2

fc[32309,2]=3  #131 counts | AluSc8SINE/Alu
fc[195667,2]=2.5  #1757 counts | L2cLINE/L2
fc[210913,2]=4  # 580 counts | L3LINE/CR1
fc[215423,2]=2  # 25 counts  | LTR17LTR/ERV1
fc[119083,2]=1.5  # 28 counts | HERV17-intLTR/ERV1
fc[216236,2]=3 # 709 counts | LTR2BLTR/ERV1
fc[330301,2]=4  # 8 counts | Tigger3aDNA/TcMar-Tigger
fc[232141,2]=2 # 110 counts | MER20DNA/hAT-Charlie
fc[138814,2]=1.5 # 257 counts | L1MA4LINE/L1
fc[305226,2]=2 # 22 counts | MLT1A0LTR/ERVL-MaLR

# built in function of polyester to simulate reads

simulate_experiment(seqpath="data/", reads_per_transcript=counts_vect, fold_changes=fc, feature="exon", gtf="data/CP068277.2_norRNA_noSimpleRepeat_rmsk.gtf", num_reps=c(3,3), readlen=76, outdir="sim_reads2", paired=FALSE)

