suppressWarnings(library(data.table))
suppressWarnings(library(UpSetR))
suppressWarnings(library(optparse))

# local vs UGER ------------------------------------------------------------------------------------------
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
  # in a local mac, the home directory is usuaully at '/Users/[username]'
  workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
  # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
  workdir <- "/xchip/beroukhimlab/"
}

# intermediate_dir <- paste0(workdir, "/youyun/nti/analysis_files/insertions")
# # dataset = 'DIPG'
# # dataset <- "TCGA"
# dataset <- "HCMI"
# if (dataset == "DIPG") {
#   # # DIPG proximal
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_051611.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_03012300')
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_050916.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_05102316')

#   # DIPG semi-proximal
#   # insertion.sv.calls <- fread(paste0(workdir, "youyun/nti/analysis_files/insertions_SVs_processed_051611.tsv"))
#   # intermediate_dir <- paste0(intermediate_dir, "/ins_align_total_0516231428")
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_051611.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_0516231427')
# } else if (dataset == "TCGA") {
#   # # TCGA  proximal
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_030900.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_03092312')
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_050918.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_0510231628')

#   # TCGA semi-proximal
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_051518.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_0516231150')
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_051518.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_0516231148')

#   # TCGA other callers
#   # dRanger-snowman consensus
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_062215.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_0622231546/')
#   # dRanger calls
#   # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_062217.tsv'))
#   # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_0622231707/')
# } else {
#   # HCMI proximal
#   # manta
#   insertion.sv.calls <- fread(paste0(workdir, "youyun/nti/analysis_files/insertions_SVs_processed_03251503.tsv"))
#   intermediate_dir <- paste0(intermediate_dir, "/ins_align_total_0501241334/")
#   print("Manta")
#   # svaba
#   # insertion.sv.calls <- fread(paste0(workdir, "youyun/nti/analysis_files/insertions_SVs_processed_10260025.tsv"))
#   # intermediate_dir <- paste0(intermediate_dir, "/ins_align_total_1026231833/")
#   # print("SvABA")
# }

call_single_end_configuration = function(intermediate_dir, SV_file, dataset, current_time){
    intermediate_job_input_dir <- paste0(intermediate_dir, "/inputs/")
    intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
    sig_threshold <- 0.05
    print(paste0("Significance threshold: ", sig_threshold))

    insertion.sv.calls = fread(SV_file)
    insertion.sv.calls.subset <- insertion.sv.calls[ins_len <= 30 & ins_len >= 6]
    print(paste0("number of insertion SVs with insertion length 6 to 30: ", nrow(insertion.sv.calls.subset)))

    # single breakend matching calling --------------------------------------------------------------------------------------------
    match.results <- data.table(t(apply(insertion.sv.calls.subset, 1, function(x) {
        # x = as.character(insertion.sv.calls.subset[2,])
        ins_seq <- x[which(colnames(insertion.sv.calls) == "ins_seq")]
        current_breakend_ID <- x[which(colnames(insertion.sv.calls) == "breakend_ID")]
        intermediate_ins_dir <- paste0(intermediate_job_output_dir, "/alignment_files/", ins_seq)

        out_og_breakend_results <- fread(paste0(intermediate_ins_dir, "/", ins_seq, "_out_og_sig_breakends.tsv"))
        in_og_breakend_results <- fread(paste0(intermediate_ins_dir, "/", ins_seq, "_in_og_sig_breakends.tsv"))
        out_rc_breakend_results <- fread(paste0(intermediate_ins_dir, "/", ins_seq, "_out_rc_sig_breakends.tsv"))
        in_rc_breakend_results <- fread(paste0(intermediate_ins_dir, "/", ins_seq, "_in_rc_sig_breakends.tsv"))

        # Checking if out breakend_ID is in the significant match breakend list in the output in the 4 possible configurations
        out.ins.match <- c(
            out_og_breakend_results[breakend_ID == current_breakend_ID]$adjusted_quantile,
            out_og_breakend_results[breakend_ID == current_breakend_ID]$j_min
        )
        in.ins.match <- c(
            in_og_breakend_results[breakend_ID == current_breakend_ID]$adjusted_quantile,
            in_og_breakend_results[breakend_ID == current_breakend_ID]$j_min
        )
        out.ins.rc.match <- c(
            out_rc_breakend_results[breakend_ID == current_breakend_ID]$adjusted_quantile,
            out_rc_breakend_results[breakend_ID == current_breakend_ID]$j_min
        )
        in.ins.rc.match <- c(
            in_rc_breakend_results[breakend_ID == current_breakend_ID]$adjusted_quantile,
            in_rc_breakend_results[breakend_ID == current_breakend_ID]$j_min
        )

        return(c(current_breakend_ID, out.ins.match, in.ins.match, out.ins.rc.match, in.ins.rc.match))
    })))

    # merge the result back into the original sv file -----------------------------------------------------------------------
    colnames(match.results) <- c(
        "breakend_ID", "outside_ins_match_quantile", "outside_ins_match_j_max_quantile",
        "inside_ins_match_quantile", "inside_ins_match_j_max_quantile",
        "outside_ins_rc_match_quantile", "outside_ins_rc_match_j_max_quantile",
        "inside_ins_rc_match_quantile", "inside_ins_rc_match_j_max_quantile"
    )
    match.results <- match.results[, lapply(.SD, as.numeric), by = breakend_ID]
    insertion.sv.calls.aligned <- data.table(merge(insertion.sv.calls.subset, match.results, by = "breakend_ID"))

    insertion.sv.calls.aligned[, c(
        "outside_ins_match_quantile_fdr", "outside_ins_rc_match_quantile_fdr",
        "inside_ins_match_quantile_fdr", "inside_ins_rc_match_quantile_fdr"
    ) := lapply(.SD, function(x) p.adjust(x, method = "fdr")),
    .SDcols = c(
        "outside_ins_match_quantile", "outside_ins_rc_match_quantile",
        "inside_ins_match_quantile", "inside_ins_rc_match_quantile"
    )
    ]

    insertion.sv.calls.aligned[
        , c("outside_ins_match", "outside_ins_rc_match", "inside_ins_match", "inside_ins_rc_match") :=
          list(
              ifelse(outside_ins_match_quantile_fdr < sig_threshold, 1, 0),
              ifelse(outside_ins_rc_match_quantile_fdr < sig_threshold, 1, 0),
              ifelse(inside_ins_match_quantile_fdr < sig_threshold, 1, 0),
              ifelse(inside_ins_rc_match_quantile_fdr < sig_threshold, 1, 0)
          )
    ]

    if (any(grepl(insertion.sv.calls.aligned$seqnames, pattern = "chr"))) {
        insertion.sv.calls.aligned$seqnames <- factor(insertion.sv.calls.aligned$seqnames, levels = paste0("chr", c(seq(1, 22, 1), "X", "Y")))
    } else {
        insertion.sv.calls.aligned$seqnames <- factor(insertion.sv.calls.aligned$seqnames, levels = c(seq(1, 22, 1), "X", "Y"))
    }
    insertion.sv.calls.aligned[, c("chr_order", "SV_config_combo") := list(ifelse(seqnames[1] <= seqnames[2], TRUE, FALSE), paste0(cnt_type, collapse = "")),
        by = .(SV_ID = gsub(":[0-2]$", "", ID), Sample)
    ]
    print(paste0("writing insertion sv alignment calls to: ", paste0(intermediate_job_output_dir, dataset, "_insertions_w_sig_alignment_", current_time, ".tsv")))
    write.table(insertion.sv.calls.aligned, 
        paste0(intermediate_job_output_dir, dataset, "_insertions_w_sig_alignment_", current_time, ".tsv"),
        sep = "\t", row.names = FALSE
    )

    # making a upset plot on the breakend matching result
    print(paste0("writing upset plot to: ", paste0(intermediate_job_output_dir, dataset, "_alignment_upset_plot_", current_time, ".pdf")))
    pdf(paste0(intermediate_job_output_dir, dataset, "_alignment_upset_plot_", current_time, ".pdf"),
        width = 8, height = 8, onefile = FALSE
    )
    upset(insertion.sv.calls.aligned[, .(outside_ins_match, outside_ins_rc_match, inside_ins_match, inside_ins_rc_match)],
        sets.bar.color = "#56B4E9", empty.intersections = TRUE
    )
    dev.off()

    return(insertion.sv.calls.aligned)
}

call_combination_configuration = function(intermediate_dir, insertion.sv.calls.aligned, dataset, current_time){
    intermediate_job_input_dir <- paste0(intermediate_dir, "/inputs/")
    intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
    # generating the dcast data for combination information -----------------------------------------------------------------------
    # insertion.sv.calls.aligned = fread(paste0(workdir, "youyun/nti/analysis_files/HCMI_insertions_w_sig_alignment_09281828.tsv"))
    # dup and del has breakend 1 vs breakend 2
    # translocations and inversions have breakend 0 and breakend 1
    insertion.sv.calls.aligned[, c("SV_ID", "breakend_num") :=
        list(
            paste0(Sample, "__", gsub(":[0-2]$", "", ID)),
            gsub(".*:", "", ID)
        )]
    # per sv pair, we sort the breakend numbers, sort them, and give them a unified order
    insertion.sv.calls.aligned[order(breakend_num), breakend_order := c("breakend1", "breakend2"), SV_ID]
    # dcast into the short format with breakend 1 and 2 identified by SV_ID
    insertion_sv_calls_aligned_paired <- dcast(
      # melt into long format
        melt(
            insertion.sv.calls.aligned,
            id.vars = c(
                "SV_ID", "Sample", "breakend_order", "SV_config_combo",
                # "ins_count", "ins_count_bin", 
                "ins_len", "both_end_pass_filter"
            )
        )[, variable := paste0(variable, "_", breakend_order)][, breakend_order := NULL],
        SV_ID + Sample + SV_config_combo + 
        # ins_count + ins_count_bin + 
        ins_len + both_end_pass_filter ~ variable,
        value.var = "value"
    )
    nrow(insertion_sv_calls_aligned_paired)

    # comprehensively calculate each possible combination of the breakend matching result ----------------------------------------
    # we have outside_ins_match, outside_ins_rc_match, inside_ins_match, inside_ins_rc_match
    # we want to calculate every single possible combination of these 4 variables
    insertion_sv_calls_aligned_paired[,
        c(
            "outside_outside", "outside_outside_rc", "outside_inside", "outside_inside_rc",
            "outside_rc_outside_rc", "outside_rc_inside", "outside_rc_inside_rc",
            "inside_inside", "inside_inside_rc", "inside_rc_inside_rc"
        ) := list(
            ifelse(outside_ins_match_breakend1 == 1 & outside_ins_match_breakend2 == 1, 1, 0),
            ifelse((outside_ins_match_breakend1 == 1 & outside_ins_rc_match_breakend2 == 1) |
                (outside_ins_match_breakend2 == 1 & outside_ins_rc_match_breakend1 == 1), 1, 0),
            ifelse((outside_ins_match_breakend1 == 1 & inside_ins_match_breakend2 == 1) |
                (outside_ins_match_breakend2 == 1 & inside_ins_match_breakend1 == 1), 1, 0),
            ifelse((outside_ins_match_breakend1 == 1 & inside_ins_rc_match_breakend2 == 1) |
                (outside_ins_match_breakend2 == 1 & inside_ins_rc_match_breakend1 == 1), 1, 0),
            ifelse(outside_ins_rc_match_breakend1 == 1 & outside_ins_rc_match_breakend2 == 1, 1, 0),
            ifelse((outside_ins_rc_match_breakend1 == 1 & inside_ins_match_breakend2 == 1) |
                (outside_ins_rc_match_breakend2 == 1 & inside_ins_match_breakend1 == 1), 1, 0),
            ifelse((outside_ins_rc_match_breakend1 == 1 & inside_ins_rc_match_breakend2 == 1) |
                (outside_ins_rc_match_breakend2 == 1 & inside_ins_rc_match_breakend1 == 1), 1, 0),
            ifelse(inside_ins_match_breakend1 == 1 & inside_ins_match_breakend2 == 1, 1, 0),
            ifelse((inside_ins_match_breakend1 == 1 & inside_ins_rc_match_breakend2 == 1) |
                (inside_ins_match_breakend2 == 1 & inside_ins_rc_match_breakend1 == 1), 1, 0),
            ifelse(inside_ins_rc_match_breakend1 == 1 & inside_ins_rc_match_breakend2 == 1, 1, 0)
        ),
        by = .(SV_ID)
    ]
    # output the combination file
    print(paste0("writing combination file to: ", paste0(intermediate_job_output_dir, dataset, "_alignment_combination_", current_time, ".tsv")))
    write.table(insertion_sv_calls_aligned_paired, paste0(intermediate_job_output_dir, dataset, "_alignment_combination_", current_time, ".tsv"),
        sep = "\t", row.names = FALSE
    )

    # high confidence combination calling ---------------------------------------------------------------------------------------
    # I want to make sure that every row, breakend1 only has one match and breakend2 also only has one match
    insertion_sv_calls_aligned_paired[
        , c("breakend1_match", "breakend2_match") :=
            list(
                sum(as.numeric(c(outside_ins_match_breakend1, outside_ins_rc_match_breakend1, inside_ins_match_breakend1, inside_ins_rc_match_breakend1))),
                sum(as.numeric(c(outside_ins_match_breakend2, outside_ins_rc_match_breakend2, inside_ins_match_breakend2, inside_ins_rc_match_breakend2)))
            ),
        by = .(SV_ID)
    ]
    table(insertion_sv_calls_aligned_paired$breakend1_match)
    table(insertion_sv_calls_aligned_paired$breakend2_match)
    high_conf_insertion_sv_calls_aligned_paired <- insertion_sv_calls_aligned_paired[
        (outside_outside == 1 | outside_outside_rc == 1 | outside_inside == 1 | outside_inside_rc == 1 |
          outside_rc_outside_rc == 1 | outside_rc_inside == 1 | outside_rc_inside_rc == 1 |
          inside_inside == 1 | inside_inside_rc == 1 | inside_rc_inside_rc == 1) &
          (breakend1_match == 1 & breakend2_match == 1),
    ]
    print(high_conf_insertion_sv_calls_aligned_paired[
        , .(
            SV_ID, Sample, cnt_type_breakend1, cnt_type_breakend2, ins_seq_breakend1, ins_seq_breakend2,
            breakend1_match, breakend2_match,
            outside_outside, outside_outside_rc, outside_inside, outside_inside_rc,
            outside_rc_outside_rc, outside_rc_inside, outside_rc_inside_rc,
            inside_inside, inside_inside_rc, inside_rc_inside_rc
        )
    ])
}


if(!interactive()) {
    option_list = list(
        make_option(c("-o", "--outputdir"),
            type = "character",
            help = "Output directory to use.", metavar = "outputdir"
        ),
        make_option(c("-n", "--datasetname"),
            type = "character", default = "hcmi_cancer_models",
            help = "Dataset to use. Current option: [HCMI, 1KG, hcmi_cancer_models]", 
            metavar = "datasetname"
        )
    )

    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)

    outputdir = opt$outputdir
    dataset = opt$datasetname

    if (is.null(outputdir)) {
        stop("Output directory is required.")
    }

    setwd(paste0(workdir, "youyun/nti/code/insertion_SVs"))
    source('align_and_config_call_helper.R')
    current_time <- format(Sys.time(), "%m%d%H%M")


    print("Finding SV file")
    SV_file = find_vcf_file(outputdir)

    print("Calling single end configuration")
    single_end_significance = call_single_end_configuration(
      outputdir, SV_file, dataset, current_time
    )

    print("Calling double end configuration")
    double_end_significance = call_combination_configuration(
      outputdir, single_end_significance, dataset, current_time
    )
}
