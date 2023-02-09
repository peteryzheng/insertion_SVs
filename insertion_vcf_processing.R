library(data.table)
library(stringr)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

# reading input ----------
# local vs UGER
if(Sys.getenv("LOGNAME") == 'youyunzheng'){
  workdir = '~/Documents/HMS/PhD/beroukhimlab/broad_mount/'
}else{
  workdir = '/xchip/beroukhimlab/'
}

# PCAWG path
# sv_files = list.files(paste0(workdir,'siyun/data/insertions/pcawg/'),
#                       pattern = 'sv.vcf$',full.names = TRUE)
# # DIPG
sv_files = list.files(paste0(workdir,'Frank/DIPG/SVABAvcfs'),
                      pattern = 'sv.vcf$',full.names = TRUE)

# rbind them into one df
colnames9 <- c("seqnames","start","ID","REF","ALT","QUAL","FILTER","INFO","Sample")

sv.calls = do.call(rbind,lapply(sv_files,
                                function(x){
                                  # read vcf file
                                  tmp = data.table::fread(cmd=paste("grep -v '^#'", x),sep='\t')
                                  if(nrow(tmp) == 0){
                                    tmp = data.table(matrix(nrow = 0,ncol = 9))
                                  }else{
                                    tmp = tmp[,c(1:8)]
                                    # assign sample name column
                                    tmp$Sample = gsub('.*/|.svaba.*','',x)
                                  }
                                  # assign the correct column names
                                  colnames(tmp) = colnames9
                                  return(tmp)
                                }
))

# FILTERS ----------
# filter for minimal qual score
sv.calls = sv.calls[as.numeric(QUAL) > 15,]
# filter for sv breakpoints with insertion sequences only
insertion.sv.calls = sv.calls[grepl('INSERTION',INFO) & seqnames %in% c(1:22,'X','Y'),]
insertion.sv.calls[,SV_ID := paste(Sample,seqnames,ID,sep = '__')]
# making sure every breakend will have its pair (both breakends are on chr 1-22)
insertion.sv.calls[,both_end_pass_filter := ifelse(.N == 2,TRUE,FALSE),,by = .(SV_pair_ID = gsub(':.*','',ID),Sample)] 
insertion.sv.calls = insertion.sv.calls[both_end_pass_filter == TRUE]

print(paste0('Total number of SV breakends (qual score <= 15 discarded): ',nrow(sv.calls)))
print(paste0('Total number of samples before insertion filter: ',length(sv_files)))
print(paste0('Total number of samples before insertion filter (samples with any SVs): ',length(unique(sv.calls$Sample))))
print(paste0('Percentage of SV breakends with MH before insertion filter: ',round(nrow(sv.calls[grepl('HOMSEQ',INFO)])/nrow(sv.calls),3)))

print(paste0('Total number of SV breakends after filtering for insertions (and chr 1-22,X,Y): ',nrow(insertion.sv.calls)))
print(paste0('Total number of samples after filter: ',length(unique(insertion.sv.calls$Sample))))
print(paste0('Percentage of SV breakends with MH after filter: ',round(nrow(insertion.sv.calls[grepl('HOMSEQ',INFO)])/nrow(insertion.sv.calls),3)))

# Extracting info ----------
# extract insertion sequences from the info column
insertion.sv.calls$ins_seq = gsub('INSERTION=','',str_extract(insertion.sv.calls$INFO,'INSERTION=[ACGT]*'))
insertion.sv.calls$ins_len = nchar(insertion.sv.calls$ins_seq)
# extracting MH
insertion.sv.calls$mh_seq = gsub('HOMSEQ=','',str_extract(insertion.sv.calls$INFO,'HOMSEQ=[ACGT]*'))
insertion.sv.calls$mh_len = nchar(insertion.sv.calls$mh_seq)
# extracting evidence type
insertion.sv.calls$evdnc = gsub('EVDNC=','',str_extract(insertion.sv.calls$INFO,'EVDNC=[A-Z\\_]*'))
# extracting MAPQ
insertion.sv.calls$mapq = as.numeric(gsub('MAPQ=','',str_extract(insertion.sv.calls$INFO,'MAPQ=[0-9]*')))
insertion.sv.calls$disc_mapq = as.numeric(gsub('DISC_MAPQ=','',str_extract(insertion.sv.calls$INFO,'DISC_MAPQ=[0-9]*')))
# extracting +/- breakpoint
insertion.sv.calls$cnt_type = unlist(lapply(grepl('^[AGCT]',insertion.sv.calls$ALT),function(x) {ifelse(x,'+','-')}))
colnames(insertion.sv.calls)

write.table(insertion.sv.calls,paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_',format(Sys.time(), "%m%d%H"),'.tsv'),
            sep = '\t',row.names = FALSE)



# run matching ----------
# match.results = t(apply(insertion.sv.calls[1,], 1, align_nearby,max.mismatch.num = 2))
# colnames(match.results) = c('outside_ins_match','inside_ins_match','outside_ins_rc_match','inside_ins_rc_match',
#                             'outside_mh_match','inside_mh_match','outside_mh_rc_match','inside_mh_rc_match')
# insertion.sv.calls = data.table(cbind(insertion.sv.calls,match.results))




x = as.character(unlist(insertion.sv.calls[ID == '2001128:1' & Sample == '10_417-5243_pair']))
intermediate_dir = paste0(workdir,'/youyun/nti/analysis_files/insertions')
# create intermediate directory to store all intermediate alignment results for RAM efficiency
intermediate_dir = paste0(intermediate_dir,'/ins_align_',format(Sys.time(), "%m%d%y%H"))
dir.create(intermediate_dir,showWarnings = TRUE)
print(paste0('Creating SV directory: ',intermediate_dir))

system.time(align_nearby(x,insertion.sv.calls,intermediate_dir,7,1))
system.time(align_nearby_mc(x,insertion.sv.calls,intermediate_dir,7,1))

# testing multithreading 
x = as.character(unlist(insertion.sv.calls[ID == '12572961:1' & Sample == 'GBM12_pair']))
intermediate_dir = paste0(workdir,'/youyun/nti/analysis_files/insertions')
system.time(align_nearby(x,insertion.sv.calls,intermediate_dir,7,1))
system.time(align_nearby_mc(x,insertion.sv.calls,intermediate_dir,7,1))

