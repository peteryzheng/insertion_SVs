suppressWarnings(library(data.table))


# =================================== batch_align_kmers.R ===================================

create_folder_structure = function(intermediate_dir){
    current_time <- gsub('.*_|/$','',intermediate_dir)
    # file paths -----------------------------------------------------------------------------------------------------------
    # create intermediate directory to store all intermediate alignment results for RAM efficiency
    intermediate_job_input_dir <- paste0(intermediate_dir, "/inputs/")
    intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
    print(paste0("Output directory here: ", intermediate_dir))
    dir.create(intermediate_dir, showWarnings = TRUE)
    dir.create(intermediate_job_input_dir, showWarnings = TRUE)
    dir.create(intermediate_job_output_dir, showWarnings = TRUE)
    dir.create(paste0(intermediate_job_output_dir,'/task_array_output'), showWarnings = TRUE)
    dir.create(paste0(intermediate_job_output_dir,'/alignment_files'), showWarnings = TRUE)
    dir.create(paste0(intermediate_job_output_dir,'/tmp_working_directory'), showWarnings = TRUE)
}

write_kmer_file = function(intermediate_dir, SV_file){
    current_time <- gsub('.*_|/$','',intermediate_dir)
    intermediate_job_input_dir <- paste0(intermediate_dir, "/inputs/")

    insertion.sv.calls <- fread(SV_file)
    insertion.sv.calls.subset <- insertion.sv.calls[ins_len <= 30 & ins_len >= 6]

    # write the kmers to a file -----------------------------------------------------------------------------------------------------------
    kmer_file = paste0(intermediate_job_input_dir, "kmers_", current_time, ".txt")
    # adding in sample to prevent the more prevalent kmers from being run first all the time
    write.table(sample(unique(insertion.sv.calls.subset$ins_seq)), kmer_file, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
    print(paste0("The kmers are here: ", kmer_file))
    return(kmer_file)
}

# =================================== batch_align_kmers_rerun.r ===================================

check_qacct_report = function(intermediate_dir, task_array_number, kmer_taskid){
    intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
    # getting jobs info ----------------------------------------------------------------------------------------------------
    qstat_output <- system(paste0("qacct -j ", task_array_number), intern = TRUE)

    # processing jobs info for plotting run time -------------------------------------------------------------------------------------------------
    task_status <- data.table(
        taskid = as.numeric(gsub("taskid| ", "", qstat_output[grep("taskid", qstat_output)])),
        wallclock = as.numeric(gsub("ru_wallclock| ", "", qstat_output[grep("ru_wallclock", qstat_output)])) / 60 / 60,
        maxvmem = as.numeric(gsub("maxvmem| |G", "", qstat_output[grep("maxvmem", qstat_output)])),
        exit_status = gsub("exit_status| ", "", qstat_output[grep("exit_status", qstat_output)])
    )
    task_status_kmer <- merge(task_status, kmer_taskid, by = "taskid")
    print(paste0("Number of kmers that finished successfully: ", nrow(task_status_kmer[exit_status == 0])))

    if (any(task_status_kmer$exit_status != 0)) {
        print("Some kmer failed")
        # UGER is missing info
        print(paste0("Number of kmers that failed (using UGER): ", nrow(task_status_kmer[exit_status != 0])))
    } else {
        print("All kmers finished successfully (using UGER)")
    }
    print(paste0("Total number of kmers we have uger information for: ", nrow(task_status_kmer)))
    write.table(task_status_kmer, paste0(
        intermediate_job_output_dir,'/', task_array_number, ".tsv"
    ), sep = "\t", quote = FALSE, row.names = FALSE)
    print(paste0('Job resource info from UGER can be found here: ', paste0(
        intermediate_job_output_dir,'/', task_array_number, ".tsv"
    )))
    # this is for when uger can't run ggplot
    # task_status_kmer = fread('/xchip/beroukhimlab//youyun/nti/analysis_files/insertions/ins_align_total_1026231833/outputs//40825613_1031231405.tsv')
    ggplot(task_status_kmer[exit_status == 0], aes(x = factor(nchar(V1)), y = wallclock)) +
        geom_violin() +
        geom_jitter(height = 0, width = 0.1) +
        scale_x_discrete(breaks = seq(5, 31, 1))
    # print('ggplot command')
    ggsave(paste0(intermediate_job_output_dir,"runtime_plot_", task_array_number, ".pdf"), width = 10, height = 10)
    print(paste0("Runtime plot is here: ", paste0(intermediate_job_output_dir,"runtime_plot_", task_array_number, ".pdf")))
}

check_task_array_output = function(intermediate_dir, task_array_number, kmer_taskid){
    current_time = format(Sys.time(), "%m%d%y%H%M")
    intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
    # actually checking if every job finished running by the text output ------------------------------------------------------
    # looking at the actual output files
    # output_files <- list.files(paste0(workdir, "/youyun/nti/code/outputs/"), pattern = paste0(task_name, ".o", task_array_number, ".*"), full.names = TRUE)
    output_files <- list.files(
        paste0(intermediate_job_output_dir, "/task_array_output/"), 
        pattern = paste0("\\.o", task_array_number, ".*"), 
        full.names = TRUE
    )
    job_status_lookup <- data.table(do.call("rbind", lapply(output_files, function(x) {
        current_output <- readLines(x)
        # all successful jobs should have the 'Done!' in the output
        success <- ifelse(any(grepl("Done!", current_output)), TRUE, FALSE)
        num_files <- ifelse(length(current_output) >= 3, length(list.files(gsub('"', "", gsub(paste0(".*/xchip/beroukhimlab/"), workdir, current_output[3])))), 0)
        task_id <- gsub(".*\\.", "", x)
        kmer <- ifelse(length(current_output) > 0, current_output[1], kmer_taskid$V1[as.numeric(task_id)])
        num_lines <- length(current_output)
        c(kmer = kmer, task_id = task_id, success = success, num_lines = num_lines, num_files = num_files, file_name = x)
    })))
    # jobs that finished, (with either 13 lines of output or last line is 'Done!') + (16 files in the output directory)
    finished_jobs <- job_status_lookup[success == TRUE & num_files == 12]
    print(paste0("Number of jobs that finished successfully: ", nrow(finished_jobs)))
    # jobs that failed, (with less than 13 lines of output or last line is not 'Done!') + (less than 16 files in the output directory)
    failed_jobs <- job_status_lookup[success == FALSE | num_files != 12]
    print(paste0("Number of jobs that failed: ", nrow(failed_jobs)))
    print(paste0("They should add up to the total number of jobs: ", nrow(job_status_lookup)))
    if (nrow(failed_jobs) == 0) {
        # stop the code
        stop("All jobs finished successfully")
    }
    # writing the failed kmer file
    writeLines(text = failed_jobs$kmer, paste0(intermediate_dir, "/inputs/failed_kmers_", current_time, ".txt"), sep = "\n", useBytes = FALSE)
    print(paste0("The failed kmers are here: ", paste0(intermediate_dir, "/inputs/failed_kmers_", current_time, ".txt")))
    return(paste0(intermediate_dir, "/inputs/failed_kmers_", current_time, ".txt"))
}

# =================================== batch align shared functions ===================================

generate_qsub_script = function(
    intermediate_dir, alignparam, kmer_file, SV_file, script_name, 
    downsample_num, seed, time, cores
){
    current_time <- gsub('.*_|/$','',intermediate_dir)
    intermediate_job_input_dir <- paste0(intermediate_dir, "/inputs/")
    intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
    # write the task array script -----------------------------------------------------------------------------------------------------------
    # ALIGNMENT PARAMETERS -- for all the penalty values, positive and negative doesnt matter
    if(alignparam == "mhe"){
        # 33n2 SCHEME THAT THE MHe project uses
        # gap_open = 7, gap_epen = 1, mismatch_pen = 1, match_pen = 3
        align_str = "-g 7 -e 1 -m 1 -t 3"
    }else if(alignparam == "default_bwa"){
        # DEFAULT VALUE OF BWA-MEM
        # https://bio-bwa.sourceforge.net/bwa.shtml
        # gap_open = 6, gap_epen = 1, mismatch_pen = 4, match_pen = 1
        align_str = '-g 6 -e 1 -m 4 -t 1'
    }
    # INCREASEING THE MATCH REWARD TO HOPEFULLY GET HIGHER TEMPLATE FIDELITY
    # gap_open = 7, gap_epen = 1, mismatch_pen = 1, match_pen = 9
    # "-g 7 -e 1 -m 1 -t 9"

    template_task_array <- c(
        "#!/bin/bash",
        paste0("#$ -l h_rt=",time,":00:00"),
        paste0("#$ -t 1-", nrow(fread(kmer_file, header = FALSE))),
        paste0("#$ -pe smp ",cores," "),
        paste0("#$ -binding linear:",cores," "),
        "#$ -l h_vmem=4G",
        paste0("#$ -o '", intermediate_job_output_dir, "/task_array_output/' "),
        paste0("#$ -e '", intermediate_job_output_dir, "/task_array_output/' "),
        paste0("#$ -N ins_homeology_", current_time),
        "",
        'export PATH="/xchip/beroukhimlab/youyun/miniconda3/bin:$PATH"',
        "",
        "# Going into core dump directory and clean up",
        paste0('cd ',intermediate_job_output_dir,'/tmp_working_directory'),
        "rm -rf core*",
        '',
        paste0("kmer=$(sed -n -e \"$SGE_TASK_ID p\"  ", kmer_file, ")"),
        "echo $kmer",
        paste0(
            "/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n rameen Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/align_nearby_utils.R  -i $kmer ",
            # parameters for how far out the breakpoint to align to
            " -w 20 -d ", SV_file, " -o ", intermediate_job_output_dir,'/alignment_files/', " ",
            # " -b 100 -d /xchip/beroukhimlab/youyun/nti/analysis_files/", SV_file, " -o ", intermediate_dir, " "
            align_str,' -n ', downsample_num, ' -s ', seed
        )
    )
    task_array_path <- paste0(intermediate_job_input_dir, script_name, "_", current_time, ".sh")
    writeLines(text = template_task_array, task_array_path, sep = "\n", useBytes = FALSE)
    print(paste0("The task array script is here: ", task_array_path))
    # system('use UGER')
    print(paste0("qsub ", task_array_path))
}

generate_terra_file = function(intermediate_dir, alignparam, kmer_file, downsample_num, seed){
    current_time <- gsub('.*_|/$','',intermediate_dir)
    intermediate_job_input_dir <- paste0(intermediate_dir, "/inputs/")
    intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
    
    kmers = fread(kmer_file, header = FALSE)
    colnames(kmers) = c('entity:sample_id')
    
    # write the terra sample table -----------------------------------------------------------------------------------------------------------
    # ALIGNMENT PARAMETERS -- for all the penalty values, positive and negative doesnt matter
    if(alignparam == "mhe"){
        # 33n2 SCHEME THAT THE MHe project uses
        # gap_open = 7, gap_epen = 1, mismatch_pen = 1, match_pen = 3
        align_str = "-g 7 -e 1 -m 1 -t 3"
        kmers$gap_open = 7
        kmers$gap_epen = 1
        kmers$mismatch_pen = 1
        kmers$match_pen = 3
    }else if(alignparam == "default_bwa"){
        # DEFAULT VALUE OF BWA-MEM
        # https://bio-bwa.sourceforge.net/bwa.shtml
        # gap_open = 6, gap_epen = 1, mismatch_pen = 4, match_pen = 1
        align_str = '-g 6 -e 1 -m 4 -t 1'
        kmers$gap_open = 6
        kmers$gap_epen = 1
        kmers$mismatch_pen = 4
        kmers$match_pen = 1
    }

    kmers$downsample_num = downsample_num
    kmers$seed = seed

    terra_file = paste0(intermediate_job_input_dir, "terra_sample_table_", current_time, ".tsv")
    write.table(kmers, terra_file, sep = "\t", quote = FALSE, row.names = FALSE)
    print(paste0("The terra sample table is here: ", terra_file))
    

}

# =================================== batch align and config calling shared functions ===================================

find_vcf_file = function(intermediate_dir){
    current_time <- gsub('.*_|/$','',intermediate_dir)

    SV_file <- list.files(
        intermediate_dir, pattern = "insertions_SVs_processed_[0-9]*.tsv", 
        full.names = TRUE
    )
    if(length(SV_file) == 0){
        stop("No insertion SV file found in the directory.")
    }else if(length(SV_file) > 1){
        stop("More than one insertion SV file found in the directory.")
    }else {
        print(paste0("SV file used is this: ", SV_file))
    }
    return(SV_file)
}
