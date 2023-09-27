library(data.table)
library(ggplot2)

# local vs UGER
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
    # in a local mac, the home directory is usuaully at '/Users/[username]'
    workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
    # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
    workdir <- "/xchip/beroukhimlab/"
}

# inputs ---------------------------------------------------------------------------------------------------------------
intermediate_dir <- paste0(workdir, "/youyun/nti/analysis_files/insertions")
# create intermediate directory to store all intermediate alignment results for RAM efficiency
################## change this to the directory you want to store the results ##################
intermediate_dir = paste0(intermediate_dir, "/ins_align_total_0918231455")
print(paste0("Output directory here: ", intermediate_dir))
intermediate_job_output_dir = paste0(intermediate_dir,'/outputs')

# previously used SV file from the task array w failed jobs
################## change this to the SV file you want to use ##################
SV_file <- "insertions_SVs_processed_09181429.tsv"
print(paste0("SV file used is this: ", workdir, "youyun/nti/analysis_files/", SV_file))

# kmer file from the task array w failed jobs
################## change this to the kmer file you want to use ##################
kmer_text <- fread(paste0(workdir, "/youyun/nti/code/outputs/kmers_0918231455.txt"), header = FALSE)
# kmer_text <- fread(paste0(workdir, "/youyun/nti/code/outputs/kmers_0918231455.txt"), header = FALSE)

# jobs number for the task array
################## change this to the task array number you want to use ##################
task_array_number <- 40346375
task_name <- "ins_homeology_0918231455"
# task_array_number <- 40346375

# getting jobs info ----------------------------------------------------------------------------------------------------
qstat_output <- system(paste0("qacct -j ", task_array_number), intern = TRUE)

# processing jobs info for plotting run time -------------------------------------------------------------------------------------------------
kmer_text$taskid <- c(1:nrow(kmer_text))
task_status <- data.table(
    taskid = as.numeric(gsub("taskid| ", "", qstat_output[grep("taskid", qstat_output)])),
    wallclock = as.numeric(gsub("ru_wallclock| ", "", qstat_output[grep("ru_wallclock", qstat_output)])) / 60 / 60,
    maxvmem = as.numeric(gsub("maxvmem| |G", "", qstat_output[grep("maxvmem", qstat_output)])),
    exit_status = gsub("exit_status| ", "", qstat_output[grep("exit_status", qstat_output)])
)
task_status_kmer <- merge(task_status, kmer_text, by = "taskid")
print(paste0("Number of kmers that finished successfully: ", nrow(task_status_kmer[exit_status == 0])))
current_time <- format(Sys.time(), "%m%d%y%H%M")
if (any(task_status_kmer$exit_status != 0)) {
    print("Some kmer failed")
    # UGER is missing info
    print(paste0("Number of kmers that failed (using UGER): ", nrow(task_status_kmer[exit_status != 0])))
}
print(paste0("total number of kmers we have uger information for: ", nrow(task_status_kmer)))

ggplot(task_status_kmer[exit_status == 0], aes(x = factor(nchar(V1)), y = wallclock)) +
    geom_violin() +
    geom_jitter(height = 0, width = 0.1) +
    scale_x_discrete(breaks = seq(5, 31, 1))
ggsave(paste0(workdir, "/youyun/nti/analysis_files/runtime_plot_", task_array_number, "_", current_time, ".pdf"), width = 10, height = 10)
print(paste0("Runtime plot is here: ", paste0(workdir, "/youyun/nti//analysis_files/runtime_plot_", task_array_number, "_", current_time, ".pdf")))

# actually checking if every job finished running by the text output ------------------------------------------------------
# looking at the actual output files
output_files <- list.files(paste0(workdir, "/youyun/nti/code/outputs/"), pattern = paste0(task_name, ".o", task_array_number, ".*"), full.names = TRUE)
# output_files <- list.files(intermediate_job_output_dir, pattern = paste0(task_name, ".o", task_array_number, ".*"), full.names = TRUE)
job_status_lookup = data.table(do.call('rbind',lapply(output_files, function(x) {
    current_output = readLines(x)
    # all successful jobs should have 13 lines of output
    success = ifelse(length(current_output) == 13,TRUE,FALSE)
    # all successful jobs should have the 'Done!' in the output
    # success = ifelse(any(current_output == 'Done!'), TRUE, FALSE)        
    num_files = ifelse(length(current_output) >= 3,length(list.files(gsub('"','',gsub(paste0('.*/xchip/beroukhimlab/'),workdir,current_output[3])))),0)
    task_id = gsub('.*\\.', "", x)
    kmer = ifelse(length(current_output) > 0, current_output[1],kmer_text$V1[as.numeric(task_id)])
    num_lines = length(current_output)
    c(kmer = kmer, task_id = task_id, success = success,num_lines = num_lines,num_files = num_files,file_name = x)
})))
# jobs that finished, (with either 13 lines of output or last line is 'Done!') + (16 files in the output directory)
finished_jobs <- job_status_lookup[success == TRUE & num_files == 16]
print(paste0("Number of jobs that finished successfully: ", nrow(finished_jobs)))
# jobs that failed, (with less than 13 lines of output or last line is not 'Done!') + (less than 16 files in the output directory)
failed_jobs <- job_status_lookup[success == FALSE | num_files != 16]
print(paste0("Number of jobs that failed: ", nrow(failed_jobs)))
print(paste0('They should add up to the total number of jobs: ', nrow(job_status_lookup)))
# writing the failed kmer file
writeLines(text = failed_jobs$kmer, paste0(workdir, "/youyun/nti/code/outputs/failed_kmers_", current_time, ".txt"), sep = "\n", useBytes = FALSE)
print(paste0("The failed kmers are here: ", paste0(workdir, "/youyun/nti/code/outputs/failed_kmers_", current_time, ".txt")))

# submitting jobs ------------------------------------------------------------------------------------------------------
template_task_array <- c(
    "#!/bin/bash",
    "#$ -l h_rt=60:00:00",
    paste0("#$ -t 1-", nrow(fread(paste0(workdir, "/youyun/nti/code/outputs/failed_kmers_", current_time, ".txt"), header = FALSE))),
    "#$ -pe smp 4 ",
    "#$ -binding linear:4 ",
    "#$ -l h_vmem=4G",
    paste0("#$ -o '",intermediate_job_output_dir,"' "),
    paste0("#$ -e '",intermediate_job_output_dir,"' "),
    paste0("#$ -N ins_homeology_", format(Sys.time(), "%m%d%y%H%M")),
    "",
    "source /broad/software/scripts/useuse",
    "use R-4.0",
    "",
    paste0("kmer=$(sed -n -e \"$SGE_TASK_ID p\"  ", workdir, "/youyun/nti/code/outputs/failed_kmers_", current_time, ".txt", ")"),
    "echo $kmer",
    paste0(
        "Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/align_nearby_utils.R  -i $kmer ",
        # parameters for how far out the breakpoint to align to
        " -w 20 -d /xchip/beroukhimlab/youyun/nti/analysis_files/", SV_file, " -o ", intermediate_dir, " ",
        # " -b 100 -d /xchip/beroukhimlab/youyun/nti/analysis_files/", SV_file, " -o ", intermediate_dir, " "
        # alignment parameters -- for all the penalty values, positive and negative doesnt matter
        # 33n2 scheme that the MHe project uses
        # gap_open = 7, gap_epen = 1, mismatch_pen = 1, match_pen = 3
        # the default values for bwa mem
        # https://bio-bwa.sourceforge.net/bwa.shtml
        # gap_open = 6, gap_epen = 1, mismatch_pen = 4, match_pen = 1
        # '-g 6 -e 1 -m 4 -t 1'
        "-g 7 -e 1 -m 1 -t 3"
    )
)
task_array_path <- paste0(workdir, "youyun/nti/code/outputs/task_array_rerun_", current_time, ".sh")
writeLines(text = template_task_array, task_array_path, sep = "\n", useBytes = FALSE)
print(paste0("The task array script is here: ", task_array_path))
# system('use UGER')
print(paste0("qsub ", task_array_path))
