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
    # print(summary(best_alignment$quantile))
    # adjusting max qunatile
    best_alignment$adjusted_quantile <- 1 - rank(best_alignment$quantile) / length(best_alignment$quantile)
    return(best_alignment)
}

score_to_p_val = function(score_mat, breakend_IDs, kmer){
    quantile_mat = data.table(apply(score_mat, 2, function(x) rank(x) / length(x)))[, breakend_ID := breakend_IDs]
    significance_df = alignment_sig(quantile_mat, kmer)
    significance_df[!grepl('background',breakend_ID)]
}

align_nearby_mc_bp <- function(ins_seq, bases, insertion.sv.calls, intermediate_dir, gap_open, gap_epen, mismatch_pen, match_pen, ncores) {
    bases <- as.numeric(bases)
    # match penalty is the original value, mismatch penalty is the negative of the input value
    mat <- nucleotideSubstitutionMatrix(match = match_pen, mismatch = mismatch_pen, type = "DNA")[c("A", "G", "C", "T", "N"), c("A", "G", "C", "T", "N")]

    # creating intermeadiate directories to store all alignment matricies in
    intermediate_ins_dir <- paste0(intermediate_dir, "/", ins_seq)
    dir.create(intermediate_ins_dir, showWarnings = TRUE)
    print(paste0("Creating ins directory for ", ins_seq, ": ", intermediate_ins_dir))

    # old implementation of alignment scoring -- requires high vmem bc insertion.sv.calls is called in each thread?
    # ins_alignment_scores <- mclapply(X = c(1:nrow(insertion.sv.calls)), FUN = function(x) {
    #     each_call <- as.character(unlist(insertion.sv.calls[x, ]))
    #     chr <- paste0("chr", each_call[which(colnames(insertion.sv.calls) == "seqnames")])
    #     start <- as.numeric(each_call[which(colnames(insertion.sv.calls) == "start")])
    #     cnt <- each_call[which(colnames(insertion.sv.calls) == "cnt_type")]
    #     out_refseq <- each_call[which(colnames(insertion.sv.calls) == "outside_ref")]
    #     in_refseq <- each_call[which(colnames(insertion.sv.calls) == "inside_ref")]
    #     # return(list(find_best_alignment(DNAString(ins_seq),window,chr,start,cnt,'outside',gap_open,gap_epen,mat),
    #     #             find_best_alignment(DNAString(ins_seq),window,chr,start,cnt,'inside',gap_open,gap_epen,mat),
    #     #             find_best_alignment(reverseComplement(DNAString(ins_seq)),window,chr,start,cnt,'outside',gap_open,gap_epen,mat),
    #     #             find_best_alignment(reverseComplement(DNAString(ins_seq)),window,chr,start,cnt,'inside',gap_open,gap_epen,mat)))
    #     return(list(
    #         find_best_alignment_substring(DNAString(ins_seq), bases, chr, cnt, "outside", out_refseq, gap_open, gap_epen, mat),
    #         find_best_alignment_substring(DNAString(ins_seq), bases, chr, cnt, "inside", in_refseq, gap_open, gap_epen, mat),
    #         find_best_alignment_substring(reverseComplement(DNAString(ins_seq)), bases, chr, cnt, "outside", out_refseq, gap_open, gap_epen, mat),
    #         find_best_alignment_substring(reverseComplement(DNAString(ins_seq)), bases, chr, cnt, "inside", in_refseq, gap_open, gap_epen, mat)
    #     ))
    # }, mc.cores = ncores, mc.preschedule = TRUE, mc.set.seed = 55555)

    # new implementation of alignment scoring -- hopefully requires less vmem
    ins_alignment_scores = mcmapply(
        function(chr,start,cnt,out_refseq,in_refseq){
            return(list(
                find_best_alignment_substring(DNAString(ins_seq), bases, chr, cnt, "outside", out_refseq, gap_open, gap_epen, mat),
                find_best_alignment_substring(DNAString(ins_seq), bases, chr, cnt, "inside", in_refseq, gap_open, gap_epen, mat),
                find_best_alignment_substring(reverseComplement(DNAString(ins_seq)), bases, chr, cnt, "outside", out_refseq, gap_open, gap_epen, mat),
                find_best_alignment_substring(reverseComplement(DNAString(ins_seq)), bases, chr, cnt, "inside", in_refseq, gap_open, gap_epen, mat)
            ))
        }, 
        insertion.sv.calls$seqnames, insertion.sv.calls$start, insertion.sv.calls$cnt_type, 
        insertion.sv.calls$outside_ref, insertion.sv.calls$inside_ref,
        mc.cores = ncores, mc.preschedule = TRUE, mc.set.seed = 55555,
        USE.NAMES = FALSE
    )
    print(dim(ins_alignment_scores))

    # Insertion alignment score matrix ----------
    ins_alignment_scores_out_og <- data.table(do.call("rbind", ins_alignment_scores[1,]))
    ins_alignment_scores_in_og <- data.table(do.call("rbind", ins_alignment_scores[2,]))
    ins_alignment_scores_out_rc <- data.table(do.call("rbind", ins_alignment_scores[3,]))
    ins_alignment_scores_in_rc <- data.table(do.call("rbind", ins_alignment_scores[4,]))

    print(paste0("Writing alignment scores to ", intermediate_ins_dir))
    write.table(ins_alignment_scores_out_og, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_out_og.tsv"), sep = "\t", row.names = FALSE)
    write.table(ins_alignment_scores_in_og, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_in_og.tsv"), sep = "\t", row.names = FALSE)
    write.table(ins_alignment_scores_out_rc, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_out_rc.tsv"), sep = "\t", row.names = FALSE)
    write.table(ins_alignment_scores_in_rc, paste0(intermediate_ins_dir, "/", ins_seq, "_alignment_score_in_rc.tsv"), sep = "\t", row.names = FALSE)

    # column wise quantile matrix ----------
    # ins_alignment_scores_out_og <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_og.tsv'),sep = '\t')
    # ins_alignment_scores_in_og <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_og.tsv'),sep = '\t')
    # ins_alignment_scores_out_rc <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_rc.tsv'),sep = '\t')
    # ins_alignment_scores_in_rc <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_rc.tsv'),sep = '\t')

    real_sv_breakend_IDs = grep('background',insertion.sv.calls$breakend_ID, value = TRUE, invert = TRUE)
    background_sv_breakend_IDs = grep('background',insertion.sv.calls$breakend_ID, value = TRUE)
    cat(paste0(
        'Out of ',nrow(insertion.sv.calls),' breakends, \n',
        '\t ',length(real_sv_breakend_IDs),' are real insertion SVs, \n',
        '\t ',length(background_sv_breakend_IDs),' are background breakends'
    ))

    print("Calculating significance values for outside og")
    background_mat_out_og = ins_alignment_scores_out_og[grepl('background',insertion.sv.calls$breakend_ID),]
    print(dim(ins_alignment_scores_out_og))
    print(dim(background_mat_out_og))
    out.ins.match = do.call('rbind',mclapply(real_sv_breakend_IDs,function(x){
        # for each real insertion SV, calculate the p value separately with the background breakends
        score_to_p_val(
            # row bind to create a [N_background + 1] x [bases out to align] matrix
            rbind(background_mat_out_og, ins_alignment_scores_out_og[which(insertion.sv.calls$breakend_ID == x)]),
            c(background_sv_breakend_IDs,x), ins_seq
        )
    }, mc.cores = ncores, mc.preschedule = FALSE, mc.set.seed = 55555))
    write.table(out.ins.match, paste0(intermediate_ins_dir, "/", ins_seq, "_out_og_sig_breakends.tsv"), sep = "\t", row.names = FALSE)
    
    print("Calculating significance values for inside og")
    background_mat_in_og = ins_alignment_scores_in_og[grepl('background',insertion.sv.calls$breakend_ID),]
    print(dim(ins_alignment_scores_in_og))
    print(dim(background_mat_in_og))
    in.ins.match = do.call('rbind',mclapply(real_sv_breakend_IDs,function(x){
        score_to_p_val(
            rbind(background_mat_in_og, ins_alignment_scores_in_og[which(insertion.sv.calls$breakend_ID == x)]),
            c(background_sv_breakend_IDs,x), ins_seq
        )
    }, mc.cores = ncores, mc.preschedule = FALSE, mc.set.seed = 55555))
    write.table(in.ins.match, paste0(intermediate_ins_dir, "/", ins_seq, "_in_og_sig_breakends.tsv"), sep = "\t", row.names = FALSE)

    print("Calculating significance values for outside rc")
    background_mat_out_rc = ins_alignment_scores_out_rc[grepl('background',insertion.sv.calls$breakend_ID),]
    print(dim(ins_alignment_scores_out_rc))
    print(dim(background_mat_out_rc))
    out.ins.rc.match = do.call('rbind',mclapply(real_sv_breakend_IDs,function(x){
        score_to_p_val(
            rbind(background_mat_out_rc, ins_alignment_scores_out_rc[which(insertion.sv.calls$breakend_ID == x)]),
            c(background_sv_breakend_IDs,x), ins_seq
        )
    }, mc.cores = ncores, mc.preschedule = FALSE, mc.set.seed = 55555))
    write.table(out.ins.rc.match, paste0(intermediate_ins_dir, "/", ins_seq, "_out_rc_sig_breakends.tsv"), sep = "\t", row.names = FALSE)

    print("Calculating significance values for inside rc")
    background_mat_in_rc = ins_alignment_scores_in_rc[grepl('background',insertion.sv.calls$breakend_ID),]
    print(dim(ins_alignment_scores_in_rc))
    print(dim(background_mat_in_rc))
    in.ins.rc.match = do.call('rbind',mclapply(real_sv_breakend_IDs,function(x){
        score_to_p_val(
            rbind(background_mat_in_rc, ins_alignment_scores_in_rc[which(insertion.sv.calls$breakend_ID == x)]),
            c(background_sv_breakend_IDs,x), ins_seq
        )
    }, mc.cores = ncores, mc.preschedule = FALSE, mc.set.seed = 55555))
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
        ),
        make_option(c("-n", "--downsamplenum"),
            type = "numeric", default = 4000,
            help = "number of breakends to downsample to", metavar = "downsample_num"
        ),
        make_option(c("-s", "--seed"),
            type = "numeric", default = 55555,
            help = "seed for random number generator", metavar = "seed"
        ),
        make_option(c("-c", "--cores"),
            type = "numeric", default = 4,
            help = "Number of cores to use", metavar = "cores"
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
    downsample_num <- as.numeric(opt$downsamplenum)
    seed <- as.numeric(opt$seed)
    ncores <- as.numeric(opt$cores)

    if(intermediate_dir == '.'){
        intermediate_dir = getwd()
    }

    set.seed(seed)
    
    # downsample the breakends
    if(nrow(insertion.sv.calls) > downsample_num){
        print(paste0(
            "Downsampling to ", downsample_num, " breakends"
        ))
        background_indices = sample(1:nrow(insertion.sv.calls), downsample_num)
        insertion.sv.calls <- insertion.sv.calls[
            c(
                # background
                background_indices,
                # real insertion SVs
                which(insertion.sv.calls$ins_seq == ins)
            ), 
        ][
            ,breakend_ID := paste0(
                breakend_ID, c(
                    # give a tag to the background breakend IDs
                    paste0('_background_',c(1:downsample_num)),
                    # the real insertion SV breakend IDs keep the same
                    rep('', sum(insertion.sv.calls$ins_seq == ins))
                )
            )
        ]
    }else{
        print(paste0(
            "Using all ", nrow(insertion.sv.calls), " breakends"
        ))
        insertion.sv.calls <- insertion.sv.calls[
            ,breakend_ID := paste0(
                breakend_ID, ifelse(
                    ins_seq == ins, '', '_background'
                )
            )
        ]
    }
    print(nrow(insertion.sv.calls))

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
    print(paste0(
        "Aligning the insertion sequence [", ins, "] with a ", total_bases, 
        "bp window using all breakends in file [", opt$data, 
        "] and writing to [", intermediate_dir, "]"
    ))
    print(paste0(
        "Alignment parameters: gap_open = ", gapopen, ", gap_epen = ", gapepen, 
        ", mismatch_pen = ", mismatchpen, ", match_pen = ", matchpen
    ))

    system.time(unlist(align_nearby_mc_bp(
        ins, total_bases, insertion.sv.calls, intermediate_dir,
        # aligmnent parameters from make_option inputs
        gap_open = gapopen, gap_epen = gapepen, mismatch_pen = mismatchpen, match_pen = matchpen,
        ncores = ncores
    )))
    print('Done!')
}