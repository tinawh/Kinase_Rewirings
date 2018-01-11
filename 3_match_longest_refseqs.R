
require(data.table); require(tidyr); require(dplyr); require(reshape2)

loc <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/"

# files
load("/.mounts/labs/reimandlab/private/users/thuang/data/05-29-17-take1/all_RefGene_proteins.fa.rsav") #fa1

aachanges <- readRDS("/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/aachanges.rds")
full.mut <- aachanges$AAChange.refGene # take
# file names 
longest.refseqs.file <- "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/longest_refseqs.rds"

# functions
ExtractInfo <- function(row) {
    nm <- unlist(regmatches(row, gregexpr("NM_[0123456789]+", row)))
    p <- unlist(regmatches(row, gregexpr("p.[abcdefghijklmnopqrstuvwxyz0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ]+", row)))
    p <- substr(p, 3, nchar(p))
    paste(nm, p)
 }

LongestRefSeq <- function(row) {
	# find length of each id from fasta file and remove all except longest ones from all_refmuts 
    row <- unlist(row)
    # single variant 
    if (length(row) == 1) {
        row
    # multiple variants
    } else {
        refs <- unname(nchar(sapply(row, function(x) fa1[unlist(regmatches(x, gregexpr("NM_[0123456789]+", x)))])))
        i <- which(refs == max(refs, na.rm = TRUE))
        if (length(i) == 1) {
            row[which.max(refs)]
        } else {   
            max(unlist(sapply(i, function(x) row[x])), na.rm = TRUE)
        }  
    }
}
    
Organize <- function(longest.refs) {
	binded <- cbind(aachanges, longest.refs, stringsAsFactors = F) 
	binded.with.names <- binded %>% 
        					mutate(hugo_name = gsub(":.*", "", AAChange.refGene)) %>% 
        					rename(patient = Otherinfo, longest_refs = longest.refs)
    return(binded.with.names)
}
####################################################################################################################################
all.refmuts <- unname(sapply(full.mut, ExtractInfo))
longest.refs <- sapply(all.refmuts, LongestRefSeq)
binded.with.names <- Organize(longest.refs)

saveRDS(binded.with.names, "/.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/longest_refseqs.rds")

