suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(parallel))

# Helper functions ----------

# given a SV breakend (location, +/-, base pair, alignment penalty parameters) and outside/inside, align insertions to k-2k bp of reference sequence and get (k+1) alignment scores
find_best_alignment_substring <- function(ins_seq, bases, chr, cnt, side, ref_seq, gap_open_pen, gap_ext_pen, sub_mat) {
  if ((cnt == "+" & side == "outside") || (cnt == "-" & side == "inside")) {
    # extension towards higher genomic position
    end <- seq(nchar(ins_seq), bases, 1)
    align_scores <- unlist(lapply(substring(ref_seq, 1, end), function(x) {
      pairwiseAlignment(ins_seq, DNAString(x),
        type = "global-local",
        substitutionMatrix = sub_mat, gapOpening = gap_open_pen, gapExtension = gap_ext_pen, scoreOnly = TRUE
      )
    }))
  } else if ((cnt == "-" & side == "outside") || (cnt == "+" & side == "inside")) {
    # extension towards lower genomic position
    end <- seq(nchar(ref_seq) - nchar(ins_seq) + 1, nchar(ref_seq) - bases + 1, -1)
    align_scores <- unlist(lapply(substring(ref_seq, end, nchar(ref_seq)), function(x) {
      pairwiseAlignment(ins_seq, DNAString(x),
        type = "global-local",
        substitutionMatrix = sub_mat, gapOpening = gap_open_pen, gapExtension = gap_ext_pen, scoreOnly = TRUE
      )
    }))
  }
  return(align_scores)
}

# getting significance value from a alignment score matrix
alignment_sig <- function(alignment_df, kmer) {
  # kmer = 'TTTTGTTTAATTTAAATTAAAC'
  # alignment_df = fread('/Users/youyunzheng/Documents/HMS/PhD/beroukhimlab/broad_mount//youyun/nti/analysis_files/insertions/ins_align_total_02102315/TTTTGTTTAATTTAAATTAAAC/TTTTGTTTAATTTAAATTAAAC_alignment_quantile_out_og.tsv')
  # get the best alignment score (by quantile) for every breakend to get the adjusted p value
  best_alignment <-
    data.table(melt(alignment_df, id.vars = "breakend_ID", variable.name = "j", value.name = "quantile"))[, k_mer := kmer][
      # setting j to the correct value (column 1 is supposed to be k by k, column 2 is supposed to be k by k+1. etc)
      , .(
        j = as.numeric(gsub("V", "", j)) + nchar(kmer) - 1, quantile, k_mer,
        # tagging max quantile value per breakend (breakend_ID), max quantile in k by j at the breakend
        max_quantile = ifelse(quantile == max(quantile), TRUE, FALSE)
      ),
      by = "breakend_ID"
    ][(max_quantile),
      # getting minimum avg and max offset with the max quantile value
      .(j_min = min(j), j_avg = mean(j), j_max = max(j)),
      by = c("breakend_ID", "k_mer", "max_quantile", "quantile")
    ]
  # best quantile per SV
  print(summary(best_alignment$quantile))
  # adjusting max qunatile
  best_alignment$adjusted_quantile <- rank(best_alignment$quantile) / length(best_alignment$quantile)
  # adding significance info back into the total df
  best_alignment$significance <- ifelse(best_alignment$adjusted_quantile > 0.95, TRUE, FALSE)
  return(best_alignment)
}

align_nearby_mc_bp <- function(ins_seq, bases, insertion.sv.calls, intermediate_dir, gap_open, gap_epen, mismatch_pen, match_pen) {
  bases <- as.numeric(bases)
  # match penalty is the original value, mismatch penalty is the negative of the input value
  mat <- nucleotideSubstitutionMatrix(match = match_pen, mismatch = mismatch_pen, type = "DNA")[c("A", "G", "C", "T", "N"), c("A", "G", "C", "T", "N")]

  # creating intermeadiate directories to store all alignment matricies in
  intermediate_ins_dir <- paste0(intermediate_dir, "/", ins_seq)
  dir.create(intermediate_ins_dir, showWarnings = TRUE)
  print(paste0("Creating ins directory for ", ins_seq, ": ", intermediate_ins_dir))
  ins_alignment_scores <- mclapply(X = c(1:nrow(insertion.sv.calls)), FUN = function(x) {
    each_call <- as.character(unlist(insertion.sv.calls[x, ]))
    chr <- paste0("chr", each_call[which(colnames(insertion.sv.calls) == "seqnames")])
    start <- as.numeric(each_call[which(colnames(insertion.sv.calls) == "start")])
    cnt <- each_call[which(colnames(insertion.sv.calls) == "cnt_type")]
    out_refseq <- each_call[which(colnames(insertion.sv.calls) == "outside_ref")]
    in_refseq <- each_call[which(colnames(insertion.sv.calls) == "inside_ref")]
    # return(list(find_best_alignment(DNAString(ins_seq),window,chr,start,cnt,'outside',gap_open,gap_epen,mat),
    #             find_best_alignment(DNAString(ins_seq),window,chr,start,cnt,'inside',gap_open,gap_epen,mat),
    #             find_best_alignment(reverseComplement(DNAString(ins_seq)),window,chr,start,cnt,'outside',gap_open,gap_epen,mat),
    #             find_best_alignment(reverseComplement(DNAString(ins_seq)),window,chr,start,cnt,'inside',gap_open,gap_epen,mat)))
    return(list(
      find_best_alignment_substring(DNAString(ins_seq), bases, chr, cnt, "outside", out_refseq, gap_open, gap_epen, mat),
      find_best_alignment_substring(DNAString(ins_seq), bases, chr, cnt, "inside", in_refseq, gap_open, gap_epen, mat),
      find_best_alignment_substring(reverseComplement(DNAString(ins_seq)), bases, chr, cnt, "outside", out_refseq, gap_open, gap_epen, mat),
      find_best_alignment_substring(reverseComplement(DNAString(ins_seq)), bases, chr, cnt, "inside", in_refseq, gap_open, gap_epen, mat)
    ))
  }, mc.cores = 4, mc.preschedule = TRUE, mc.set.seed = 55555)

  # Insertion alignment score matrix ----------
  ins_alignment_scores_out_og <- do.call("rbind", lapply(ins_alignment_scores, function(x) x[[1]]))
  ins_alignment_scores_in_og <- do.call("rbind", lapply(ins_alignment_scores, function(x) x[[2]]))
  ins_alignment_scores_out_rc <- do.call("rbind", lapply(ins_alignment_scores, function(x) x[[3]]))
  ins_alignment_scores_in_rc <- do.call("rbind", lapply(ins_alignment_scores, function(x) x[[4]]))

  write.table(ins_alignment_scores_out_og, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_out_og.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_scores_in_og, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_in_og.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_scores_out_rc, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_out_rc.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_scores_in_rc, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_in_rc.tsv"), sep = "\t", row.names = FALSE)

  # column wise quantile matrix ----------
  # ins_alignment_scores_out_og <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_og.tsv'),sep = '\t')
  # ins_alignment_scores_in_og <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_og.tsv'),sep = '\t')
  # ins_alignment_scores_out_rc <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_rc.tsv'),sep = '\t')
  # ins_alignment_scores_in_rc <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_rc.tsv'),sep = '\t')

  ins_alignment_quantile_out_og <- data.table(apply(ins_alignment_scores_out_og, 2, function(x) rank(x) / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]
  ins_alignment_quantile_in_og <- data.table(apply(ins_alignment_scores_in_og, 2, function(x) rank(x) / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]
  ins_alignment_quantile_out_rc <- data.table(apply(ins_alignment_scores_out_rc, 2, function(x) rank(x) / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]
  ins_alignment_quantile_in_rc <- data.table(apply(ins_alignment_scores_in_rc, 2, function(x) rank(x) / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]

  ins_alignment_quantile_out_og_max <- data.table(apply(ins_alignment_scores_out_og, 2, function(x) rank(x, ties.method = "max") / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]
  ins_alignment_quantile_in_og_max <- data.table(apply(ins_alignment_scores_in_og, 2, function(x) rank(x, ties.method = "max") / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]
  ins_alignment_quantile_out_rc_max <- data.table(apply(ins_alignment_scores_out_rc, 2, function(x) rank(x, ties.method = "max") / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]
  ins_alignment_quantile_in_rc_max <- data.table(apply(ins_alignment_scores_in_rc, 2, function(x) rank(x, ties.method = "max") / length(x)))[, breakend_ID := insertion.sv.calls$breakend_ID]

  write.table(ins_alignment_quantile_out_og, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_out_og.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_quantile_in_og, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_in_og.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_quantile_out_rc, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_out_rc.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_quantile_in_rc, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_in_rc.tsv"), sep = "\t", row.names = FALSE)

  write.table(ins_alignment_quantile_out_og_max, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_maxtie_out_og.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_quantile_in_og_max, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_maxtie_in_og.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_quantile_out_rc_max, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_maxtie_out_rc.tsv"), sep = "\t", row.names = FALSE)
  write.table(ins_alignment_quantile_in_rc_max, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_quantile_maxtie_in_rc.tsv"), sep = "\t", row.names = FALSE)

  # row max quantile value calculation ----------
  # ins_alignment_quantile_out_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_og.tsv'),sep = '\t')
  # ins_alignment_quantile_in_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_og.tsv'),sep = '\t')
  # ins_alignment_quantile_out_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_rc.tsv'),sep = '\t')
  # ins_alignment_quantile_in_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_rc.tsv'),sep = '\t')

  out.ins.match <- alignment_sig(ins_alignment_quantile_out_og, ins_seq)
  in.ins.match <- alignment_sig(ins_alignment_quantile_in_og, ins_seq)
  out.ins.rc.match <- alignment_sig(ins_alignment_quantile_out_rc, ins_seq)
  in.ins.rc.match <- alignment_sig(ins_alignment_quantile_in_rc, ins_seq)

  write.table(out.ins.match, paste0(intermediate_ins_dir, "/", ins_seq, "_out_og_sig_breakends.tsv"), sep = "\t", row.names = FALSE)
  write.table(in.ins.match, paste0(intermediate_ins_dir, "/", ins_seq, "_in_og_sig_breakends.tsv"), sep = "\t", row.names = FALSE)
  write.table(out.ins.rc.match, paste0(intermediate_ins_dir, "/", ins_seq, "_out_rc_sig_breakends.tsv"), sep = "\t", row.names = FALSE)
  write.table(in.ins.rc.match, paste0(intermediate_ins_dir, "/", ins_seq, "_in_rc_sig_breakends.tsv"), sep = "\t", row.names = FALSE)
}


if (!interactive()) {
  option_list <- list(
    make_option(c("-i", "--ins"),
      type = "character", default = "",
      help = "insertion sequence of interest to align to background", metavar = "insertion"
    ),
    make_option(c("-w", "--window"),
      type = "numeric", default = "",
      help = "value for window to search around the breakend (window size = N x [insertion length])", metavar = "window"
    ),
    # make_option(c("-b", "--bases"),
    #   type = "numeric", default = "",
    #   help = "value for base pairs to search around the breakend", metavar = "bases"
    # ),
    make_option(c("-d", "--data"),
      type = "character", default = "/xchip/beroukhimlab/youyun/nti/analysis_files/insertions_SVs_processed_020716.tsv",
      help = "file path to the total insertion SV file (processed using insertion_vcf_processing.R)", metavar = "data"
    ),
    make_option(c("-o", "--output"),
      type = "character", default = "",
      help = "file path to write intermediate results to", metavar = "output"
    ),
    # aligner paramters
    # for all the penalty values, positive and negative doesnt matter
    # 33n2 scheme that the MHe project uses
    # gap_open = 7, gap_epen = 1, mismatch_pen = 1, match_pen = 3
    # the default values for bwa mem
    # https://bio-bwa.sourceforge.net/bwa.shtml
    # gap_open = 6, gap_epen = 1, mismatch_pen = 4, match_pen = 1
    make_option(c("-g", "--gapopen"),
      type = "numeric", default = 7,
      help = "gap open penalty for alignment", metavar = "gap_open"
    ),
    make_option(c("-e", "--gapepen"),
      type = "numeric", default = 1,
      help = "gap extension penalty for alignment", metavar = "gap_epen"
    ),
    make_option(c("-m", "--mismatchpen"),
      type = "numeric", default = 1,
      help = "mismatch penalty for alignment", metavar = "mismatch_pen"
    ),
    make_option(c("-t", "--matchpen"),
      type = "numeric", default = 3,
      help = "match penalty for alignment", metavar = "match_pen"
    )
  )

  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  ins <- opt$ins
  window <- as.numeric(opt$window)
  # bases <- as.numeric(opt$bases)
  insertion.sv.calls <- fread(opt$data)
  intermediate_dir <- opt$output
  # aligner paramters
  gapopen <- as.numeric(opt$gapopen)
  gapepen <- as.numeric(opt$gapepen)
  mismatchpen <- as.numeric(opt$mismatchpen)
  matchpen <- as.numeric(opt$matchpen)

  # error checking for inputs
  if (!is.numeric(window) || is.null(ins) || nrow(insertion.sv.calls) == 0 || intermediate_dir == "") {
    print("check inputs!!")
  }
  # error checking for alignment parameters
  if (!is.numeric(gapopen) || !is.numeric(gapepen) || !is.numeric(mismatchpen) || !is.numeric(matchpen)) {
    print("check alignment parameters!!")
  }

  # using bases instead of windows
  # total_bases = bases
  total_bases = window * nchar(ins)
  print(paste0("Aligning the insertion sequence [", ins, "] with a ", total_bases, "bp window using all breakends in file [", opt$data, "] and writing to [", intermediate_dir, "]"))
  system.time(unlist(align_nearby_mc_bp(ins, total_bases, insertion.sv.calls, intermediate_dir,
      # aligmnent parameters from make_option inputs
      gap_open = gapopen, gap_epen = gapepen, mismatch_pen = mismatchpen, match_pen = matchpen
  )))
  print('Done!')
}