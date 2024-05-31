suppressWarnings(library(optparse))


# local vs UGER
if (Sys.getenv("HOME") %in%  c('/Users/youyun','/Users/youyunzheng')) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}

if(!interactive()) {
    option_list = list(
        make_option(c("-d", "--dataset"),
            type = "character", default = "hcmi_cancer_models",
            help = "Dataset to use. Current option: [HCMI, 1KG, hcmi_cancer_models]", metavar = "dataset"
        ),
        make_option(c("-c", "--caller"),
            type = "character", default = "manta",
            help = "Caller to use. Current option: [manta, svaba]", metavar = "caller"
        ),
        make_option(c("-p", "--dedup"),
            type = "logical", default = FALSE,
            help = "Whether to dedup the SVs.", metavar = "dedup"
        ),
        make_option(c("-o", "--outputdir"),
            type = "character", default = paste0(workdir, "youyun/nti/analysis_files/insertions"),
            help = "Output directory to use.", metavar = "outputdir"
        )
    )
    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)

    dataset = opt$dataset
    caller = opt$caller
    dedup = opt$dedup
    outputdir = opt$outputdir
    current_time <- format(Sys.time(), "%m%d%y%H%M")
    outputdir <- paste0(outputdir, "/ins_align_total_", current_time)
    print(paste0("Using dataset [", dataset, "] and caller [", caller, "]"))

    setwd(paste0(workdir, "youyun/nti/code/insertion_SVs"))
    source('helper_functions.R')
    suppressPackageStartupMessages(library(data.table))
    suppressPackageStartupMessages(library(ggplot2))
    suppressPackageStartupMessages(library(RColorBrewer))
    suppressPackageStartupMessages(library(stringr))
    suppressPackageStartupMessages(library(BSgenome))
    suppressPackageStartupMessages(library(ggrepel))

    # DATA ------------------
    sv_files = locate_data_files(dataset, caller)
    sv_calls = read_vcf_from_sv_files(sv_files)
    # load different references for different datasets
    load_reference_genome(dataset)

    # FILTERS AND REFORMATTING ----------
    # filter for major chromosomes
    # filter for minimal qual score (svaba), FILTER == 'PASS' (manta)
    # expand and reformat 1 row SVs (mantaDEL/DUP/INV) into 2 rows per SV 
    # make sure that every SV have both ends passing the above criteria
    sv_calls = filter_and_reformat_sv_calls(sv_calls, caller)

    # Extracting info ----------
    sv_calls = extract_info(sv_calls)

    # Deduping events for germline samples ----------
    if (dedup) {
        print('Deduping SVs')
        sv_calls = dedup_calls(sv_calls)
        print(paste0("After deduping, the number of SV breakends is: ", nrow(sv_calls)))
    }

    # put in surrounding sequence for each insertion sv breakend for downstream analysis ----------
    insertion_sv_calls <- sv_calls[!is.na(ins_seq)]
    insertion_sv_calls = add_surrounding_seq(insertion_sv_calls, 600)

    # some cohort stats ----------
    print('COHORT STATS: ')
    print(paste0("Total number of SV breakends after filtering: ", nrow(sv_calls)))
    print(paste0(
        "Total number of insertion SV breakends after filtering: ", nrow(insertion_sv_calls),
        " (", round(nrow(insertion_sv_calls) / nrow(sv_calls), 3), ")"
    ))
    print(paste0("Total number of samples: ", length(sv_files)))
    print(paste0("Total number of samples after filtering (samples with any SVs): ", length(unique(sv_calls$Sample))))
    print(paste0("Total number of samples after filtering with insertion SVs: ", length(unique(insertion_sv_calls$Sample))))
    print(paste0("Percentage of SV breakends with MH after filtering: ", round(nrow(sv_calls[!is.na(mh_len) | mh_len > 0]) / nrow(sv_calls), 3)))
    print(paste0("Percentage of insertion SV breakends with MH after filtering: ", round(nrow(insertion_sv_calls[!is.na(mh_len) | mh_len > 0]) / nrow(insertion_sv_calls), 3)))

    # OUTPUTS --------------------
    print(paste0("Output directory here: ", outputdir))
    dir.create(outputdir, showWarnings = FALSE)
    write_outputs(insertion_sv_calls, sv_calls, current_time, outputdir)

    # some visualization ----------
    visualizing_kmers(insertion_sv_calls, current_time, outputdir)
    ins_length_distribution(insertion_sv_calls, sv_calls, current_time, outputdir)
}