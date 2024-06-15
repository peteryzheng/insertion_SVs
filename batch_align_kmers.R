suppressWarnings(library(data.table))
suppressWarnings(library(optparse))

# local vs UGER -----------------------------------------------------------------------------------------------------------
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}

if(!interactive()) {
    option_list = list(
        make_option(c("-a", "--alignparam"),
            type = "character", default = "mhe",
            help = "Alignment Parameters to Use. Current option: [mhe, default_bwa]", metavar = "alignparam"
        ),
        make_option(c("-o", "--outputdir"),
            type = "character",
            help = "Output directory to use.", metavar = "outputdir"
        ),
        make_option(c("-n", "--downsamplenum"),
            type = "numeric", default = 4000,
            help = "number of breakends to downsample to", metavar = "downsample_num"
        ),
        make_option(c("-s", "--seed"),
            type = "numeric", default = 55555,
            help = "seed for random number generator", metavar = "seed"
        ),
        make_option(c('-t','--time'),
            type = 'numeric', default = 36,
            help = 'Time (hours) to run each task', metavar = 'time'
        ),
        make_option(c("-c", "--cores"),
            type = "numeric", default = 4,
            help = "Number of cores to use", metavar = "cores"
        )

    )

    opt_parser = OptionParser(option_list = option_list)
    opt = parse_args(opt_parser)

    alignparam = opt$alignparam
    outputdir = opt$outputdir
    downsample_num = opt$downsamplenum
    seed = opt$seed
    time = opt$time
    cores = opt$cores

    setwd(paste0(workdir, "youyun/nti/code/insertion_SVs"))
    source('helper_align_and_config.R')

    if (! alignparam %in% c("mhe", "default_bwa")) {
        stop("Invalid alignment parameter option.")
    }

    if (is.null(outputdir)) {
        stop("Output directory is required.")
    }

    create_folder_structure(outputdir)
    SV_file = find_vcf_file(outputdir)
    kmer_file = write_kmer_file(outputdir, SV_file)
    generate_qsub_script(
        outputdir, alignparam, kmer_file, SV_file, 
        'task_array', downsample_num, seed, time, cores
    )
    generate_terra_file(
        outputdir, alignparam, kmer_file, 
        downsample_num, seed
    )
}