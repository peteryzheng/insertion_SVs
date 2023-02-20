library(data.table)

# local vs UGER
if(Sys.getenv("LOGNAME") == 'youyunzheng'){
  workdir = '~/Documents/HMS/PhD/beroukhimlab/broad_mount/'
}else{
  workdir = '/xchip/beroukhimlab/'
}

intermediate_dir = paste0(workdir,'/youyun/nti/analysis_files/insertions')
# create intermediate directory to store all intermediate alignment results for RAM efficiency
intermediate_dir = paste0(intermediate_dir,'/ins_align_total_',format(Sys.time(), "%m%d%y%H"))
dir.create(intermediate_dir,showWarnings = TRUE)

# DIPG
insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_020716.tsv'))
# PCAWG
# insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_021614.tsv'))

status_report = data.table(do.call('rbind',lapply(unique(insertion.sv.calls[ins_len <= 30 & ins_len >= 6]$ins_seq),function(x){
  if(nchar(x) > 18){
    time = '72:00:00'
  }else if(nchar(x) > 11){
    time = '24:00:00'
  }else{
    time = '12:00:00'
  }
  qsub_command = paste0(
    # skeleton of the command -- 
    # 'Rscript /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/align_nearby_utils.R  -i TACACATATA -w 2 -d /xchip/beroukhimlab/youyun/nti/analysis_files/insertions_SVs_processed_020716.tsv -o /xchip/beroukhimlab/youyun/nti/analysis_files/insertions/ins_align_02072316'
    "qsub -l h_rt=",time," -pe smp 4 -binding linear:4  -l h_vmem=4G -o /xchip/beroukhimlab/youyun/nti/code/outputs -j y /xchip/beroukhimlab/youyun/util/basic_qsub_r.sh R-4.0 ",
    " \" /xchip/beroukhimlab/youyun/nti/code/insertion_SVs/align_nearby_utils.R  -i ",x,
    " -w 2 -d /xchip/beroukhimlab/youyun/nti/analysis_files/insertions_SVs_processed_020716.tsv -o ",intermediate_dir,"\" "
  )
  c(output = system(qsub_command,intern = TRUE),kmer = x,command = qsub_command)
})))

print('Jobs that did not submit successfully: ')
print(status_report[!grepl('has been submitted',output),])

status_report_path = paste0('/xchip/beroukhimlab/youyun/nti/code/outputs/submission_status_',format(Sys.time(), "%m%d%y%H%M"),'.tsv')
write.table(status_report,status_report_path,sep = '\t',row.names = FALSE)
print(paste0('Check the submission fails and results here: ',status_report_path))

