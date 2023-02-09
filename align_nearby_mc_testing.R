library(data.table)
library(stringr)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(parallel)

# Helper functions ----------
# given a SV breakend (location, +/-, window, alignment penalty parameters) and outside/inside, align insertions to k-2k bp of reference sequence and get (k+1) alignment scores
find_best_alignment = function(ins_seq,window_extend,chr,start,cnt,side,gap_open,gap_epen){
  if((cnt == '+' & side == 'outside') || (cnt == '-' & side == 'inside')){
    # extension towards higher genomic position
    end = start+nchar(ins_seq)+c(0:window_extend)
    align_scores = unlist(lapply(end,function(x) {
      pairwiseAlignment(ins_seq,getSeq(Hsapiens,chr,start,x),type = 'global-local',scoreOnly = TRUE)
    }))
  }else if((cnt == '-' & side == 'outside') || (cnt == '+' & side == 'inside')){
    # extension towards lower genomic position
    end = start-nchar(ins_seq)-c(0:window_extend)
    align_scores = unlist(lapply(end,function(x) {
      pairwiseAlignment(ins_seq,getSeq(Hsapiens,chr,x,start),type = 'global-local',scoreOnly = TRUE)
    }))
  }
  return(align_scores)
}

alignment_sig = function(alignment_df){
  # get the best alignment score (by quantile) for every breakend to get the adjusted p value
  best.alignments = apply(alignment_df[,-c('SV_ID')],1,max)
  print(summary(best.alignments))
  adjusted.best.alignments = rank(best.alignments)/length(best.alignments)
  return(alignment_df[which(adjusted.best.alignments > 0.95),]$SV_ID)
}


# align nearby methods ---------
# given a SV call, generate alignment matricies for the SV call and their respective background
align_nearby = function(x,insertion.sv.calls,intermediate_dir,gap_open,gap_epen){
  # create intermediate directory to store all intermediate alignment results for RAM efficiency
  intermediate_dir = paste0(intermediate_dir,'/ins_align_',format(Sys.time(), "%m%d%y%H"))
  dir.create(intermediate_dir,showWarnings = TRUE)
  print(paste0('Creating SV directory: ',intermediate_dir))
  ins_seq = x[12]
  mh_seq = x[14]
  chr = paste0('chr',x[1])
  start = as.numeric(x[2])
  cnt = x[19]
  SV_ID = x[10]
  window = nchar(ins_seq)
  out.ins.match <- in.ins.match <- out.ins.rc.match <- in.ins.rc.match <- out.mh.match <- in.mh.match <- out.mh.rc.match <- in.mh.rc.match <- 0
  
  # creating intermeadiate directories to store all alignment matricies in
  intermediate_ins_dir = paste0(intermediate_dir,'/',ins_seq)
  dir.create(intermediate_ins_dir,showWarnings = TRUE)
  print(paste0('Creating ins directory for ',ins_seq,': ',intermediate_ins_dir))
  insertion.sv.calls_test = insertion.sv.calls[1:10,]
  # ins_alignment_scores = apply(insertion.sv.calls,1,function(each_call) {
  ins_alignment_scores = apply(insertion.sv.calls_test,1,function(each_call) {
    list(find_best_alignment(DNAString(ins_seq),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'outside',gap_open,gap_epen),
         find_best_alignment(DNAString(ins_seq),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'inside',gap_open,gap_epen),
         find_best_alignment(reverseComplement(DNAString(ins_seq)),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'outside',gap_open,gap_epen),
         find_best_alignment(reverseComplement(DNAString(ins_seq)),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'inside',gap_open,gap_epen))
    
  })
  ins_alignment_scores_out_og = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[1]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  ins_alignment_scores_in_og = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[2]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  ins_alignment_scores_out_rc = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[3]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  ins_alignment_scores_in_rc = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[4]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  
  write.table(ins_alignment_scores_out_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_out_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_rc.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_rc.tsv'),sep = '\t',row.names = FALSE)
  
  # ins_alignment_scores_out_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_og.tsv'),sep = '\t')
  # ins_alignment_scores_in_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_og.tsv'))
  # ins_alignment_scores_out_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_rc.tsv'))
  # ins_alignment_scores_in_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_rc.tsv'))
  # 
  # out.ins.match = alignment_sig(ins_alignment_scores_out_og)
  # in.ins.match = alignment_sig(ins_alignment_scores_in_og)
  # out.ins.rc.match = alignment_sig(ins_alignment_scores_out_rc)
  # in.ins.rc.match = alignment_sig(ins_alignment_scores_in_rc)
  # c(out.ins.match,in.ins.match,out.ins.rc.match,in.ins.rc.match)
  # 
  # return(c(out.ins.match,in.ins.match,out.ins.rc.match,in.ins.rc.match,out.mh.match,in.mh.match,out.mh.rc.match,in.mh.rc.match))
}

# given a SV call, generate alignment matricies for the SV call and their respective background
align_nearby_mc = function(x,insertion.sv.calls,intermediate_dir,gap_open,gap_epen){
  # create intermediate directory to store all intermediate alignment results for RAM efficiency
  intermediate_dir = paste0(intermediate_dir,'/ins_align_',format(Sys.time(), "%m%d%y%H"))
  dir.create(intermediate_dir,showWarnings = TRUE)
  print(paste0('Creating SV directory: ',intermediate_dir))
  ins_seq = x[12]
  mh_seq = x[14]
  chr = paste0('chr',x[1])
  start = as.numeric(x[2])
  cnt = x[19]
  SV_ID = x[10]
  window = nchar(ins_seq)
  out.ins.match <- in.ins.match <- out.ins.rc.match <- in.ins.rc.match <- out.mh.match <- in.mh.match <- out.mh.rc.match <- in.mh.rc.match <- 0
  
  # creating intermeadiate directories to store all alignment matricies in
  intermediate_ins_dir = paste0(intermediate_dir,'/',ins_seq,'_mc')
  dir.create(intermediate_ins_dir,showWarnings = TRUE)
  print(paste0('Creating ins directory for ',ins_seq,': ',intermediate_ins_dir))
  insertion.sv.calls_test = insertion.sv.calls[1:10,]
  # ins_alignment_scores = apply(X = c(1:nrow(insertion.sv.calls)),1,function(each_call) {
  ins_alignment_scores = mclapply(X = c(1:nrow(insertion.sv.calls_test)),FUN = function(x) {
    each_call = as.character(unlist(insertion.sv.calls[x,]))
    return(list(find_best_alignment(DNAString(ins_seq),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'outside',gap_open,gap_epen),
                find_best_alignment(DNAString(ins_seq),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'inside',gap_open,gap_epen),
                find_best_alignment(reverseComplement(DNAString(ins_seq)),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'outside',gap_open,gap_epen),
                find_best_alignment(reverseComplement(DNAString(ins_seq)),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'inside',gap_open,gap_epen)))
  },mc.cores = 8,mc.preschedule = TRUE,mc.set.seed = 55555)
  ins_alignment_scores_out_og = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[1]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  ins_alignment_scores_in_og = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[2]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  ins_alignment_scores_out_rc = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[3]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  ins_alignment_scores_in_rc = data.table(apply(do.call('rbind',lapply(ins_alignment_scores,function(x) x[[4]])),2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls_test$SV_ID]
  
  write.table(ins_alignment_scores_out_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_out_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_rc.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_rc.tsv'),sep = '\t',row.names = FALSE)

  # ins_alignment_scores_out_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_og.tsv'),sep = '\t')
  # ins_alignment_scores_in_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_og.tsv'))
  # ins_alignment_scores_out_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_out_rc.tsv'))
  # ins_alignment_scores_in_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_in_rc.tsv'))
  # 
  # out.ins.match = alignment_sig(ins_alignment_scores_out_og)
  # in.ins.match = alignment_sig(ins_alignment_scores_in_og)
  # out.ins.rc.match = alignment_sig(ins_alignment_scores_out_rc)
  # in.ins.rc.match = alignment_sig(ins_alignment_scores_in_rc)
  # c(out.ins.match,in.ins.match,out.ins.rc.match,in.ins.rc.match)
  # 
  # return(c(out.ins.match,in.ins.match,out.ins.rc.match,in.ins.rc.match,out.mh.match,in.mh.match,out.mh.rc.match,in.mh.rc.match))

}

# print()


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

# run matching ----------
# match.results = t(apply(insertion.sv.calls[1,], 1, align_nearby,max.mismatch.num = 2))
# colnames(match.results) = c('outside_ins_match','inside_ins_match','outside_ins_rc_match','inside_ins_rc_match',
#                             'outside_mh_match','inside_mh_match','outside_mh_rc_match','inside_mh_rc_match')
# insertion.sv.calls = data.table(cbind(insertion.sv.calls,match.results))



detectCores()
x = as.character(unlist(insertion.sv.calls[ID == '2001128:1' & Sample == '10_417-5243_pair']))
intermediate_dir = paste0(workdir,'/youyun/nti/analysis_files/insertions')
system.time(align_nearby(x,insertion.sv.calls,intermediate_dir,7,1))
system.time(align_nearby_mc(x,insertion.sv.calls,intermediate_dir,7,1))

