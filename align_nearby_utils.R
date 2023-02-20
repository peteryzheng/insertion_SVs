suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(parallel))

# Helper functions ----------
# given a SV breakend (location, +/-, window, alignment penalty parameters) and outside/inside, align insertions to k-2k bp of reference sequence and get (k+1) alignment scores
find_best_alignment = function(ins_seq,window_extend,chr,start,cnt,side,gap_open,gap_epen){
  if((cnt == '+' & side == 'outside') || (cnt == '-' & side == 'inside')){
    # extension towards higher genomic position
    end = start+nchar(ins_seq)+c(0:window_extend)-1
    align_scores = unlist(lapply(end,function(x) {
      pairwiseAlignment(ins_seq,getSeq(Hsapiens,chr,start,x),type = 'global-local',scoreOnly = TRUE)
    }))
  }else if((cnt == '-' & side == 'outside') || (cnt == '+' & side == 'inside')){
    # extension towards lower genomic position
    end = start-nchar(ins_seq)-c(0:window_extend)+1
    align_scores = unlist(lapply(end,function(x) {
      pairwiseAlignment(ins_seq,getSeq(Hsapiens,chr,x,start),type = 'global-local',scoreOnly = TRUE)
    }))
  }
  return(align_scores)
}

# getting significance value from a alignment score matrix
alignment_sig = function(alignment_df){
  # alignment_df = fread(alignment_path)
  # get the best alignment score (by quantile) for every breakend to get the adjusted p value
  best_alignment = data.table(melt(alignment_df,id.vars = 'SV_ID',variable.name = 'j',value.name = 'quantile'))[
    # tagging max quantile value per sv ID
    ,.(j,quantile,max_quantile = ifelse(quantile == max(quantile),TRUE,FALSE)),by = 'SV_ID'
  ]
  # best quantile per SV
  best_alignment_per_SV = unique(best_alignment[(max_quantile),.(SV_ID,quantile)])
  print(summary(best_alignment_per_SV$quantile))
  # adjusting max qunatile
  best_alignment_per_SV$adjusted_quantile = rank(best_alignment_per_SV$quantile)/length(best_alignment_per_SV$quantile)
  # adding significance info back into the total df
  best_alignment$significance = ifelse(best_alignment$SV_ID %in% best_alignment_per_SV[adjusted_quantile > 0.95]$SV_ID,TRUE,FALSE)
  return(best_alignment)
}

# align nearby function with multicore process
# given a SV call, generate alignment matricies for the SV call and their respective background
align_nearby_mc = function(ins_seq,window,insertion.sv.calls,intermediate_dir,gap_open,gap_epen){
  window = nchar(ins_seq)*as.numeric(window)
  # out.ins.match <- in.ins.match <- out.ins.rc.match <- in.ins.rc.match <- out.mh.match <- in.mh.match <- out.mh.rc.match <- in.mh.rc.match <- 0
  
  # creating intermeadiate directories to store all alignment matricies in
  intermediate_ins_dir = paste0(intermediate_dir,'/',ins_seq)
  dir.create(intermediate_ins_dir,showWarnings = TRUE)
  print(paste0('Creating ins directory for ',ins_seq,': ',intermediate_ins_dir))
  ins_alignment_scores = mclapply(X = c(1:nrow(insertion.sv.calls)),FUN = function(x) {
    each_call = as.character(unlist(insertion.sv.calls[x,]))
    return(list(find_best_alignment(DNAString(ins_seq),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'outside',gap_open,gap_epen),
                find_best_alignment(DNAString(ins_seq),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'inside',gap_open,gap_epen),
                find_best_alignment(reverseComplement(DNAString(ins_seq)),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'outside',gap_open,gap_epen),
                find_best_alignment(reverseComplement(DNAString(ins_seq)),window,paste0('chr',each_call[1]),as.numeric(each_call[2]),each_call[19],'inside',gap_open,gap_epen)))
  },mc.cores = 4,mc.preschedule = TRUE,mc.set.seed = 55555)

  ins_alignment_scores_out_og = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[1]]))
  ins_alignment_scores_in_og = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[2]]))
  ins_alignment_scores_out_rc = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[3]]))
  ins_alignment_scores_in_rc = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[4]]))

  ins_alignment_quantile_out_og = data.table(apply(ins_alignment_scores_out_og,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_in_og = data.table(apply(ins_alignment_scores_in_og,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_out_rc = data.table(apply(ins_alignment_scores_out_rc,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_in_rc = data.table(apply(ins_alignment_scores_in_rc,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]

  write.table(ins_alignment_scores_out_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_out_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_rc.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_rc.tsv'),sep = '\t',row.names = FALSE)

  write.table(ins_alignment_quantile_out_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_in_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_out_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_rc.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_in_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_rc.tsv'),sep = '\t',row.names = FALSE)

  ins_alignment_quantile_out_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_og.tsv'),sep = '\t')
  ins_alignment_quantile_in_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_og.tsv'),sep = '\t')
  ins_alignment_quantile_out_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_rc.tsv'),sep = '\t')
  ins_alignment_quantile_in_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_rc.tsv'),sep = '\t')

  out.ins.match = alignment_sig(ins_alignment_quantile_out_og)
  in.ins.match = alignment_sig(ins_alignment_quantile_in_og)
  out.ins.rc.match = alignment_sig(ins_alignment_quantile_out_rc)
  in.ins.rc.match = alignment_sig(ins_alignment_quantile_in_rc)

  write.table(out.ins.match,paste0(intermediate_ins_dir,'/',ins_seq,'_out_og_sig_breakends.tsv'),sep = '\t',row.names = FALSE)
  write.table(in.ins.match,paste0(intermediate_ins_dir,'/',ins_seq,'_in_og_sig_breakends.tsv'),sep = '\t',row.names = FALSE)
  write.table(out.ins.rc.match,paste0(intermediate_ins_dir,'/',ins_seq,'_out_rc_sig_breakends.tsv'),sep = '\t',row.names = FALSE)
  write.table(in.ins.rc.match,paste0(intermediate_ins_dir,'/',ins_seq,'_in_rc_sig_breakends.tsv'),sep = '\t',row.names = FALSE)
  
}

if(!interactive()){
  
  option_list = list(
    make_option(c("-i", "--ins"), type="character", default='', 
                help="insertion sequence of interest to align to background", metavar="insertion"),
    make_option(c("-w", "--window"), type="numeric", default="", 
                help="value for window to search around the breakend (window size = N x [insertion length])", metavar="window"),
    make_option(c("-d", "--data"), type="character", default="/xchip/beroukhimlab/youyun/nti/analysis_files/insertions_SVs_processed_020716.tsv", 
                help="file path to the total insertion SV file (processed using insertion_vcf_processing.R)", metavar="data"),
    make_option(c("-o", "--output"), type="character", default="", 
                help="file path to write intermediate results to", metavar="output")
  ); 
  
  opt_parser = OptionParser(option_list=option_list);
  opt = parse_args(opt_parser);
  
  ins = opt$ins
  window = opt$window
  insertion.sv.calls = fread(opt$data)
  intermediate_dir = opt$output
  
  if(!is.numeric(window) || is.null(ins) || nrow(insertion.sv.calls) == 0 || intermediate_dir == ''){
    print('check inputs!!')
  }
  
  # print(head(insertion.sv.calls))
  
  print(paste0('Aligning the insertion sequence [',ins,'] with a ',window,'X window using all breakends in file [',opt$data,'] and writing to [',intermediate_dir,']'))
  system.time(unlist(align_nearby_mc(ins,window,insertion.sv.calls,intermediate_dir,7,1)))
}