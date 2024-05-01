library(data.table)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(BSgenome)
library(ggrepel)

# local vs UGER
if (Sys.getenv("HOME") %in%  c('/Users/youyun','/Users/youyunzheng')) {
  # in a local mac, the home directory is usuaully at '/Users/[username]'
  workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
  # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
  workdir <- "/xchip/beroukhimlab/"
}

# Helper functions ----------

# Given a SV breakend, extract reference sequence of a certain window size inside and outside
find_surrounding_seq <- function(bases_2_extend, chr, start, cnt, side) {
  if ((cnt == "+" & side == "outside") || (cnt == "-" & side == "inside")) {
    # extension towards higher genomic position
    end <- start + bases_2_extend - 1
    if (cnt == "+" & side == "outside") {
      # outside sequence starts one bp off of the breakend
      # inside is 'inclusive' of the breakend bp
      start <- start + 1
      end <- end + 1
    }
    refseq <- getSeq(Hsapiens, chr, start, end)
  } else if ((cnt == "-" & side == "outside") || (cnt == "+" & side == "inside")) {
    # extension towards lower genomic position
    end <- start - bases_2_extend + 1
    if (cnt == "-" & side == "outside") {
      # outside sequence starts one bp off of the breakend
      # inside is 'inclusive' of the breakend bp
      start <- start - 1
      end <- end - 1
    }
    refseq <- getSeq(Hsapiens, chr, end, start)
  }
  return(as.character(refseq))
}

# DATA ----------

# PCAWG -------------------------------------
# dataset <- "pcawg"
# PCAWG SvABA
# sv_files <- list.files(paste0(workdir, "/siyun/data/insertions/pcawg/"),
#   pattern = "sv.vcf$", full.names = TRUE
# )
# consensus <- "svaba"
# PCAWG dranger
# sv_files <- list.files(paste0(workdir, "/youyun/nti/data/tcga/dr_vcf"),
#   pattern = "vcf$", full.names = TRUE
# )
# caller <- "dranger"
# PCAWG delly
# sv_files = list.files(paste0(workdir,'/youyun/nti/data/tcga/delly_vcf'),
#                       pattern = 'vcf$',full.names = TRUE)
# caller = 'delly'

# PCAWG consensus
# PCAWG total consensus
# sv_files = list.files(paste0(workdir,'/youyun/nti/data/tcga/vcf'),
#                       pattern = 'vcf$',full.names = TRUE)
# PCAWG dranger and snowman consensus
# sv_files = list.files(paste0(workdir,'/youyun/nti/data/tcga/dr_sm_vcf'),
#                       pattern = 'vcf$',full.names = TRUE)
# caller = 'consensus'

# Looking at the actual breakdown of PCAWG data in terms of caller specific calls
# sv_files = list.files(paste0(workdir,'/Simona/SV_Homeology/Data/pdcData/actualData/'),
#                       full.names = TRUE)
# sv_files = list.files(paste0(workdir,'/Simona/SV_Homeology/Data/pdcData/actualData/'),
#                       pattern = 'svfix2_4',full.names = TRUE)
# delly_files = list.files(paste0(workdir,'/Simona/SV_Homeology/Data/pdcData/actualData/'),
#                       pattern = 'delly.*.somatic',full.names = TRUE)
# dr_sm_files = list.files(paste0(workdir,'/Simona/SV_Homeology/Data/pdcData/actualData/'),
#                       pattern = 'broad-dRanger_snowman.*.somatic',full.names = TRUE)
# dr_files = list.files(paste0(workdir,'/Simona/SV_Homeology/Data/pdcData/actualData/'),
#                       pattern = 'broad-dRanger.[0-9]*.somatic',full.names = TRUE)
# sm_files = list.files(paste0(workdir,'/Simona/SV_Homeology/Data/pdcData/actualData/'),
#                       pattern = 'broad-snowman.[0-9]*.somatic',full.names = TRUE)
# print(paste0('Total number of samples: ',length(sv_files)))
# print(paste0('Total number of delly files: ',length(delly_files)))
# print(paste0('Total number of dranger snowman files: ',length(dr_sm_files)))
# print(paste0('Total number of dranger files: ',length(dr_files)))
# print(paste0('Total number of snowman files: ',length(sm_files)))

# DIPG -------------------------------------
# dataset <- "dipg"
# sv_files = list.files(paste0(workdir,'Frank/DIPG/SVABAvcfs'),
#                      pattern = 'sv.vcf$',full.names = TRUE)
# consensus = 'svaba'

# HCMI -------------------------------------
dataset <- "hcmi"
# manta
sv_files <- list.files(paste0(workdir, "youyun/nti/data/HCMI/manta"),
  pattern = "vcf$", full.names = TRUE
)
caller <- "manta"
# svaba
# sv_files <- list.files(paste0(workdir, "youyun/nti/data/HCMI/svaba"),
#   pattern = "vcf$", full.names = TRUE
# )
# caller <- "svaba"

dataset <- "hcmi_cancer_models"
# manta
sv_files <- list.files(paste0(workdir, "youyun/nti/data/HCMI/manta"),
  pattern = "-85..final.SV.WGS.vcf$", full.names = TRUE
)
caller <- "manta"

# reading input ----------
# rbind them into one df
colnames9 <- c("seqnames", "start", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "Sample")

# read in all the vcf files
sv.calls <- do.call(rbind, lapply(
  sv_files,
  function(x) {
    # read vcf file
    tmp <- data.table::fread(cmd = paste("grep -v '^#'", x), sep = "\t")
    if (nrow(tmp) == 0) {
      tmp <- data.table(matrix(nrow = 0, ncol = 9))
    } else {
      tmp <- tmp[, c(1:8)]
      # assign sample name column
      tmp$Sample <- gsub(".*/|.svaba.*|.svfix.*|.broad-dRanger_snowman.*|.final.SV.WGS.vcf", "", x)
    }
    # assign the correct column names
    colnames(tmp) <- colnames9
    return(tmp)
  }
))


# FILTERS AND REFORMATTING ----------

print(paste0("Using dataset [", dataset, "] and caller [", caller, "]"))
print(paste0('Total number of breakends: ',nrow(sv.calls)))
# load different references for different datasets
# filter for minimal qual score, samples with tumor in normal contamination
# reformat certain VCFs
if (dataset == "dipg") {
  sv.calls <- sv.calls[!Sample %in% c(
    "DIPG19_TOR_pair", "DIPG58_TOR_pair",
    "SJHGG008_A_STJ_WGS_pair", "SJHGG041_D_STJ_WGS_pair"
  ), ]
  library(BSgenome.Hsapiens.UCSC.hg19)
} else if (dataset == "pcawg") {
  library(BSgenome.Hsapiens.UCSC.hg19)
} else if (grepl("hcmi",dataset)) {
  # only HCMI uses hg38
  library(BSgenome.Hsapiens.UCSC.hg38)
}
print(head(sort(table(sv.calls[
  !seqnames %in% c(c(1:22, "X", "Y"), paste0("chr", c(1:22, "X", "Y"))) &
  !grepl('decoy|random', seqnames), 
]$seqnames), decreasing = TRUE)))
# filter for only the major chromosomes
sv.calls <- sv.calls[seqnames %in% c(c(1:22, "X", "Y"), paste0("chr", c(1:22, "X", "Y"))), ]
print(paste0('Total number of breakends after filtering for chr 1-22,X,Y,M: ',nrow(sv.calls)))
# filter for only chromosomes in the hg38 assembly in R
# sv.calls = sv.calls[
#   seqnames %in% c(names(Hsapiens@single_sequences),gsub('^chr','',names(Hsapiens@single_sequences)))
# ]
# print(paste0('Total number of breakends after filtering for all chromosomes in hg38: ',nrow(sv.calls)))

if (caller %in% c("consensus")) {
  sv.calls$callers <- gsub("CALLER=", "", str_extract(sv.calls$INFO, "CALLER=[SD]*"))
  print(table(sv.calls$callers))
  print(paste0("taking out DS or SD calls which makes up ", round(sum(sv.calls$callers == "DS" | sv.calls$callers == "SD") / nrow(sv.calls), 3), " of the data"))
  sv.calls <- sv.calls[grepl("CALLER=SD|CALLER=DS", INFO)]
} else if (caller %in% c("manta")) {
  # processing manta deletion and duplication calls
  # use ID column instead of alt column bc some deletion calls, the alt column doesn't have the mantaDEL tag
  # manta output contains deletions > 50bp that are not annotated with <DEL> in the ALT column 
  # so ID column is best 
  # sv.calls[grepl('DEL',ID) & ALT != '<DEL>',]
  dup_del_sv <- sv.calls[grepl("TANDEM|DEL", ID)]
  dup_del_sv_2nd <- copy(dup_del_sv)
  dup_del_sv[, c("ALT", "ID") := list(ifelse(
    grepl("MantaDEL", ID),
    # for the first breakend, for deletions, the lower genomic breakpoint is the + end
    # therefore they should be given t[p[ notation
    paste0(REF, "[", seqnames, ":", gsub("END=", "", str_extract(INFO, "END=[0-9]*")), "["),
    paste0("]", seqnames, ":", gsub("END=", "", str_extract(INFO, "END=[0-9]*")), "]", REF)
  ), paste0(ID, ":1"))]
  # extract the second breakend genomic coordinate for the same set of calls
  dup_del_sv_2nd[, REF := as.character(getSeq(
    Hsapiens, seqnames, as.numeric(gsub("END=", "", str_extract(INFO, "END=[0-9]*"))),
    as.numeric(gsub("END=", "", str_extract(INFO, "END=[0-9]*")))
  ))][, c("ALT", "ID") := list(ifelse(
    grepl("MantaDEL", ID),
    # for the second breakend, for deletions, the upper genomic breakpoint is the - end
    # therefore they should be given t[p[ notation
    paste0("]", seqnames, ":", start, "]", REF),
    paste0(REF, "[", seqnames, ":", start, "[")
  ), paste0(ID, ":2"))]
  dup_del_sv_2nd$start <- as.numeric(gsub("END=", "", str_extract(dup_del_sv$INFO, "END=[0-9]*")))
  sv.calls <- rbind(sv.calls[!grepl("TANDEM|DEL", ALT)], dup_del_sv, dup_del_sv_2nd)
} else {
  # consensus and manta data has no qual scores
  sv.calls <- sv.calls[as.numeric(QUAL) > 15]
}

sv.calls[, breakend_ID := paste(Sample, seqnames, ID, sep = "__")]
# making sure every breakend will have its pair (both breakends are on chr 1-22)
sv.calls[, both_end_pass_filter := ifelse(.N == 2, TRUE, FALSE), , by = .(SV_pair_ID = gsub(":[0-2]$", "", ID), Sample)]
sv.calls[, breakend_num_per_sv_ID := .N, , by = .(SV_pair_ID = gsub(":[0-2]$", "", ID), Sample)]
print(paste0("Checking if both ends pass filter: ", nrow(sv.calls)))
print(paste0("Both end pass filter: ", nrow(sv.calls[both_end_pass_filter == TRUE])))
sv.calls <- sv.calls[both_end_pass_filter == TRUE]

# Extracting info ----------

# extract insertion sequences from the info column
sv.calls$ins_seq <- gsub("INSERTION=|FORSEQ=|SVINSSEQ=", "", str_extract(sv.calls$INFO, "INSERTION=[ACGT]*|FORSEQ=[ACGT]*|SVINSSEQ=[ACGT]*"))
sv.calls$ins_len <- nchar(sv.calls$ins_seq)
# extracting MH
sv.calls$mh_seq <- gsub("HOMSEQ=", "", str_extract(sv.calls$INFO, "HOMSEQ=[ACGT]*"))
sv.calls$mh_len <- nchar(sv.calls$mh_seq)
# extracting evidence type
sv.calls$evdnc <- gsub("EVDNC=", "", str_extract(sv.calls$INFO, "EVDNC=[A-Z\\_]*"))
# extracting MAPQ
sv.calls$mapq <- as.numeric(gsub("MAPQ=", "", str_extract(sv.calls$INFO, "MAPQ=[0-9]*")))
sv.calls$disc_mapq <- as.numeric(gsub("DISC_MAPQ=", "", str_extract(sv.calls$INFO, "DISC_MAPQ=[0-9]*")))
# extracting +/- breakpoint
sv.calls$cnt_type <- unlist(lapply(grepl("^[AGCT]", sv.calls$ALT), function(x) {
  ifelse(x, "+", "-")
}))

# filter for sv breakpoints with insertion sequences only (if one breakend has insertion in info, the paired breakend will also have insertion in info) # nolint
insertion.sv.calls <- sv.calls[grepl("INSERTION|FORSEQ|SVINSSEQ", INFO), ]

print(paste0("Total number of SV breakends (qual score <= 15 discarded when qual score is available): ", nrow(sv.calls)))
# -4 bc 4 samples have tumor in normal contamination
print(paste0("Total number of samples before insertion filter: ", length(sv_files)))
print(paste0("Total number of samples before insertion filter (samples with any SVs): ", length(unique(sv.calls$Sample))))
print(paste0("Percentage of SV breakends with MH before insertion filter: ", round(nrow(sv.calls[grepl("HOMSEQ", INFO)]) / nrow(sv.calls), 3)))
print(paste0(
  "Total number of SV breakends after filtering for insertions (and chr 1-22,X,Y): ", nrow(insertion.sv.calls),
  " (", round(nrow(insertion.sv.calls) / nrow(sv.calls), 3), ")"
))
print(paste0("Total number of samples after filter: ", length(unique(insertion.sv.calls$Sample))))
print(paste0("Percentage of SV breakends with MH after filter: ", round(nrow(insertion.sv.calls[grepl("HOMSEQ", INFO)]) / nrow(insertion.sv.calls), 3)))


if (any(grepl("chr", insertion.sv.calls$seqnames))) {
  insertion.sv.calls[, c("outside_ref", "inside_ref") := list(
    find_surrounding_seq(600, seqnames, start, cnt_type, "outside"),
    find_surrounding_seq(600, seqnames, start, cnt_type, "inside")
  ), by = .(breakend_ID)]
} else {
  insertion.sv.calls[, c("outside_ref", "inside_ref") := list(
    find_surrounding_seq(600, paste0("chr", seqnames), start, cnt_type, "outside"),
    find_surrounding_seq(600, paste0("chr", seqnames), start, cnt_type, "inside")
  ), by = .(breakend_ID)]
}

# OUTPUTS --------------------

current_time <- format(Sys.time(), "%m%d%H%M")
insertion.sv.calls[, ins_count := .N, by = Sample]

# WRITING OUT FILES FOR NEXT STEPS
sv.file.path <- paste0(workdir, "youyun/nti/analysis_files/total_SVs_processed_", current_time, ".tsv")
ins.file.path <- paste0(workdir, "youyun/nti/analysis_files/insertions_SVs_processed_", current_time, ".tsv")
ins.file.filter.hypermut.path <- paste0(
  workdir, "youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_",
  current_time, ".tsv"
)
write.table(sv.calls, sv.file.path, sep = "\t", row.names = FALSE)
write.table(insertion.sv.calls, ins.file.path, sep = "\t", row.names = FALSE)
write.table(insertion.sv.calls[Sample %in% sv.calls[, .N, Sample][N <= 500, ]$Sample],
  ins.file.filter.hypermut.path,
  sep = "\t", row.names = FALSE
)
print(paste0(
  "Total number of samples after hypermutation filter: ",
  length(unique(insertion.sv.calls[Sample %in% sv.calls[, .N, Sample][N <= 500, ]$Sample]$Sample))
))
print(paste0("All SV aggregated here: ", sv.file.path))
print(paste0("Just SVs with insertions aggregated here: ", ins.file.path))
print(paste0("Just SVs with insertions without high SV burden samples here: ", ins.file.filter.hypermut.path))

# LOOK AT THE DISTRIBUTION FOR AGCT
# first getting the reference genome's ACGT distribution
for (k in c(1:3)) {
  agct_ref <- data.table(do.call("rbind", lapply(paste0("chr", c(1:22, "X", "Y")), function(x) {
    data.frame(t(oligonucleotideFrequency(BSgenome.Hsapiens.UCSC.hg38[[x]], width = k)))
  })))[, lapply(.SD, sum, na.rm = TRUE)]
  # normalize the ACGT count columns by the total count of the row
  total_bases_ref <- sum(agct_ref, na.rm = TRUE)
  agct_ref[, colnames(agct_ref) := lapply(.SD, function(x) x / sum(agct_ref, na.rm = TRUE))]
  agct_ref$source <- "ref"

  # then getting the ACGT distribution of the insertion sequences
  agct_ins <- data.table(do.call("rbind", lapply(insertion.sv.calls$ins_seq, function(x) {
    data.frame(t(oligonucleotideFrequency(DNAString(x), width = k)))
  })))[, lapply(.SD, sum, na.rm = TRUE)]
  # normalize the ACGT count columns by the total count of the row
  total_bases_ins <- sum(agct_ins, na.rm = TRUE)
  agct_ins[, colnames(agct_ins) := lapply(.SD, function(x) x / sum(agct_ins, na.rm = TRUE))]
  agct_ins$source <- "ins"

  agct_ins_sample <- data.table(do.call("rbind", lapply(unique(insertion.sv.calls$Sample), function(y) {
    oligo_freq <- data.table(do.call("rbind", lapply(insertion.sv.calls[Sample == y]$ins_seq, function(x) {
      data.frame(t(oligonucleotideFrequency(DNAString(x), width = k)))
    })))[, lapply(.SD, sum, na.rm = TRUE)]
    oligo_freq <- oligo_freq / sum(oligo_freq, na.rm = TRUE)
  })))

  kmer_p_val <- p.adjust(unlist(lapply(colnames(agct_ins_sample), function(x) {
    # print(paste0("mean of samples for base(s) ", x, ": ", mean(agct_ins_sample[[x]], na.rm = TRUE)))
    # print(paste0("reference proportion of base(s) ", x, ": ", agct_ref[[x]]))
    # print(t.test(agct_ins_sample[[x]], mu = agct_ref[[x]]))
    p_val <- t.test(agct_ins_sample[[x]], mu = agct_ref[[x]])$p.value
    # if (p_val < 0.05) {
    #   print(paste0("p-value of base(s) ", x, ": ", p_val))
    # }
    return(p_val)
  })), method = "fdr")
  print(paste0(
    "The base(s) that are significantly different between insertion sequences and reference genome after fdr is/are: ",
    paste0(colnames(agct_ins_sample)[kmer_p_val < 0.05], collapse = ", "),
    " which is ", length(colnames(agct_ins_sample)[kmer_p_val < 0.05]), " out of ", length(colnames(agct_ins_sample)), " bases"
  ))

  # plot the distribution of ACGT between reference genome and insertion sequences
  # making the bars stack on top of each other
  # make the geom text stack but dodge each other
  mycolors <- colorRampPalette(brewer.pal(8, "Set2"))(ncol(agct_ref))
  ggplot(
    merge(melt(rbind(agct_ref, agct_ins), id.vars = "source"), data.table(variable = colnames(agct_ins_sample), p_val = kmer_p_val), all.x = TRUE),
    aes(x = source, y = value, fill = variable)
  ) +
    theme(legend.position = "top") +
    geom_bar(position = "stack", stat = "identity") +
    geom_text_repel(aes(label = ifelse(p_val < 0.05, paste0(variable, ": ", round(value, 3)), "")),
      position = "stack", vjust = 0, size = 4, color = "#363434", box.padding = 0.5
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10)) +
    scale_fill_manual(values = mycolors) +
    labs(x = "Source", y = "Count") +
    ggtitle(paste0("Histogram of ACGT by Source in ", dataset, " called by ", caller))

  pdf_acgt_output_path <- paste0(workdir, "youyun/nti/analysis_files/ACGT_content_", k, "mer_", current_time, ".pdf")
  print(paste0("Saving ACGT content plot to ", pdf_acgt_output_path))
  ggsave(pdf_acgt_output_path, plot = last_plot(), width = 10, height = 10, device = "pdf")
}

# bin the samples by how many insertions they have, and make different facets for each bin
# I want to make 9 bins of equal sample size from the smaller to largest number of insertions any given sample can have in our cohort
quantiles <- quantile(unique(insertion.sv.calls[, .(ins_count), by = .(Sample)])$ins_count, seq(0, 1, 0.1))
insertion.sv.calls[, ins_count_bin := cut(ins_count, breaks = quantiles, include.lowest = TRUE)]
ggplot(insertion.sv.calls, aes(x = log10(ins_len))) +
  geom_density(aes(fill = Sample), color = "gray", alpha = 0.05) +
  labs(x = "Log10 Insertion Length", y = "Density") +
  ggtitle(paste0("Insertion Lengths in ", dataset, " Called by ", caller)) +
  theme(legend.position = "none") +
  coord_cartesian(xlim = c(0, 2)) +
  facet_wrap(~ins_count_bin, scales = "free", ncol = 3)
pdf_ins_dens_bin_output_path <- paste0(workdir, "youyun/nti/analysis_files/ins_len_density_per_ins_burden_bin_", current_time, ".pdf")
print(paste0("Saving insertion length density plot by insertion burden bin to ", pdf_ins_dens_bin_output_path))
ggsave(pdf_ins_dens_bin_output_path, plot = last_plot(), width = 10, height = 10, device = "pdf")

# plot the cumulative density and probability density distribution of insertion lengths
ggplot(insertion.sv.calls, aes(x = ins_len)) +
  geom_histogram(aes(y = ..count.. / nrow(insertion.sv.calls[ins_len == 1])), binwidth = 1, fill = "#12c3d3", alpha = 0.8) +
  scale_y_continuous(
    sec.axis = sec_axis(~ . * nrow(insertion.sv.calls[ins_len == 1]), name = "Counts")
  ) +
  geom_line(stat = "ecdf", color = "#f84c18bb") +
  labs(x = "Insertion Length", y = "CDF") +
  ggtitle(paste0("Insertion Lengths in ", dataset, " Called by ", caller)) +
  coord_cartesian(xlim = c(0, 75)) +
  geom_vline(xintercept = quantile(insertion.sv.calls$ins_len, 0.75), linetype = "dashed", color = "#d31313") +
  theme(text = element_text(size = 24))
# theme_set(theme_bw(base_size = 20))
pdf_ins_dens_output_path <- paste0(workdir, "youyun/nti/analysis_files/ins_len_density_", current_time, ".pdf")
print(paste0("Saving insertion length density plot to ", pdf_ins_dens_output_path))
ggsave(pdf_ins_dens_output_path, plot = last_plot(), width = 10, height = 10, device = "pdf")

ggplot(
  melt(merge(sv.calls[, .(total_SV = .N / 2), Sample], insertion.sv.calls[, .(ins_SV = .N / 2), Sample]), id.vars = "Sample"),
  aes(x = Sample, y = value, fill = variable)
) +
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 5)) +
  scale_fill_brewer(palette = "Set1") +
  labs(x = "Sample", y = "Count") +
  ggtitle("Histogram of total_SV and ins_SV by Sample") +
  theme_set(theme_bw(base_size = 20))
pdf_output_path <- paste0(workdir, "youyun/nti/analysis_files/SV_burden_", current_time, ".pdf")
print(paste0("Saving SV burden plot to ", pdf_output_path))
ggsave(pdf_output_path, plot = last_plot(), device = "pdf")
