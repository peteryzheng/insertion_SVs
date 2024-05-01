library(data.table)

# local vs UGER -----------------------------------------------------------------------------------------------------------
if (Sys.getenv("HOME") %in% c("/Users/youyun", "/Users/youyunzheng")) {
  # in a local mac, the home directory is usuaully at '/Users/[username]'
  workdir <- "~/Documents/HMS/PhD/beroukhimlab/broad_mount/"
} else {
  # in dipg or uger, the home directory is usuaully at '/home/unix/[username]'
  workdir <- "/xchip/beroukhimlab/"
}

# file paths -----------------------------------------------------------------------------------------------------------
current_time <- format(Sys.time(), "%m%d%y%H%M")
intermediate_dir <- paste0(workdir, "/youyun/nti/analysis_files/insertions")
# create intermediate directory to store all intermediate alignment results for RAM efficiency
intermediate_dir <- paste0(intermediate_dir, "/ins_align_total_", current_time)
intermediate_job_input_dir <- paste0(intermediate_dir, "/inputs/")
intermediate_job_output_dir <- paste0(intermediate_dir, "/outputs/")
print(paste0("Output directory here: ", intermediate_dir))
dir.create(intermediate_dir, showWarnings = TRUE)
dir.create(intermediate_job_input_dir, showWarnings = TRUE)
dir.create(intermediate_job_output_dir, showWarnings = TRUE)
dir.create(paste0(intermediate_job_output_dir,'/task_array_output'), showWarnings = TRUE)
dir.create(paste0(intermediate_job_output_dir,'/alignment_files'), showWarnings = TRUE)
dir.create(paste0(intermediate_job_output_dir,'/tmp_working_directory'), showWarnings = TRUE)

# load the insertion SV file -----------------------------------------------------------------------------------------------------------
# DIPG
# SV_file <- "insertions_SVs_processed_062217.tsv"
# SV_file = 'insertions_SVs_processed_filter_hypermut_051611.tsv'
# PCAWG -- /xchip/beroukhimlab/youyun/nti/analysis_files/insertions_SVs_processed_030900.tsv
# SV_file = 'insertions_SVs_processed_051518.tsv'
# SV_file = 'insertions_SVs_processed_filter_hypermut_051518.tsv'
# HCMI
# manta
SV_file <- "insertions_SVs_processed_03251503.tsv"
# svaba
# SV_file <- "insertions_SVs_processed_10260025.tsv"
print(paste0("SV file used is this: ", workdir, "youyun/nti/analysis_files/", SV_file))
insertion.sv.calls <- fread(paste0(workdir, "youyun/nti/analysis_files/", SV_file))
insertion.sv.calls.subset <- insertion.sv.calls[ins_len <= 30 & ins_len >= 6]

commands_text_path = paste0(intermediate_job_input_dir, "kmers_", current_time, ".txt")
write.table(unique(insertion.sv.calls.subset$ins_seq), commands_text_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
print(paste0("The kmers are here: ", commands_text_path))

template_task_array <- c(
  "#!/bin/bash",
  "#$ -l h_rt=36:00:00",
  paste0("#$ -t 1-", length(unique(insertion.sv.calls.subset$ins_seq))),
  "#$ -pe smp 4 ",
  "#$ -binding linear:4 ",
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
  paste0("kmer=$(sed -n -e \"$SGE_TASK_ID p\"  ", commands_text_path, ")"),
  "echo $kmer",
  paste0(
    "/xchip/beroukhimlab/youyun/miniconda3/bin/conda run -n rameen Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/align_nearby_utils.R  -i $kmer ",
    # parameters for how far out the breakpoint to align to
    " -w 20 -d /xchip/beroukhimlab/youyun/nti/analysis_files/", SV_file, " -o ", intermediate_job_output_dir,'/alignment_files/', " ",
    # " -b 100 -d /xchip/beroukhimlab/youyun/nti/analysis_files/", SV_file, " -o ", intermediate_dir, " "
    # ALIGNMENT PARAMETERS -- for all the penalty values, positive and negative doesnt matter
    # 33n2 SCHEME THAT THE MHe project uses
    # gap_open = 7, gap_epen = 1, mismatch_pen = 1, match_pen = 3
    "-g 7 -e 1 -m 1 -t 3"
    # DEFAULT VALUE OF BWA-MEM
    # https://bio-bwa.sourceforge.net/bwa.shtml
    # gap_open = 6, gap_epen = 1, mismatch_pen = 4, match_pen = 1
    # '-g 6 -e 1 -m 2 -t 1'
    # INCREASEING THE MATCH REWARD TO HOPEFULLY GET HIGHER TEMPLATE FIDELITY
    # gap_open = 7, gap_epen = 1, mismatch_pen = 1, match_pen = 9
    # "-g 7 -e 1 -m 1 -t 9"
  )
)
task_array_path <- paste0(intermediate_job_input_dir, "task_array_", current_time, ".sh")
writeLines(text = template_task_array, task_array_path, sep = "\n", useBytes = FALSE)
print(paste0("The task array script is here: ", task_array_path))
# system('use UGER')
print(paste0("qsub ", task_array_path))
