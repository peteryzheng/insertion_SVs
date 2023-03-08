suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
suppressPackageStartupMessages(library(parallel))

# Helper functions ----------

# Given a SV breakend, extract reference sequence of a certain window size inside and outside
find_surrounding_seq = function(bases_2_extend,chr,start,cnt,side){
  if((cnt == '+' & side == 'outside') || (cnt == '-' & side == 'inside')){
    # extension towards higher genomic position
    end = start+bases_2_extend-1
    if(cnt == '+' & side == 'outside'){
      # outside sequence starts one bp off of the breakend
      # inside is 'inclusive' of the breakend bp
      start = start + 1
      end = end + 1
    }
    refseq = getSeq(Hsapiens,chr,start,end)
  }else if((cnt == '-' & side == 'outside') || (cnt == '+' & side == 'inside')){
    # extension towards lower genomic position
    end = start-bases_2_extend+1
    if(cnt == '-' & side == 'outside'){
      # outside sequence starts one bp off of the breakend
      # inside is 'inclusive' of the breakend bp
      start = start - 1
      end = end - 1
    }
    refseq = getSeq(Hsapiens,chr,end,start)
  }
  return(as.character(refseq))
}

# given a SV breakend (location, +/-, window, alignment penalty parameters) and outside/inside, align insertions to k-2k bp of reference sequence and get (k+1) alignment scores
find_best_alignment = function(ins_seq,window_extend,chr,start,cnt,side,gap_open_pen,gap_ext_pen,sub_mat){
  if((cnt == '+' & side == 'outside') || (cnt == '-' & side == 'inside')){
    # extension towards higher genomic position
    end = start+seq(nchar(ins_seq)*1,nchar(ins_seq)*window_extend,1)-1
    if(cnt == '+' & side == 'outside'){
      # outside sequence starts one bp off of the breakend
      # inside is 'inclusive' of the breakend bp
      start = start + 1
      end = end + 1
    }
    align_scores = unlist(lapply(end,function(x) {
      refseq = getSeq(Hsapiens,chr,start,x)
      pairwiseAlignment(ins_seq,refseq,type = 'global-local',
                        substitutionMatrix = sub_mat,gapOpening = gap_open_pen,gapExtension = gap_ext_pen,scoreOnly = TRUE)
    }))
  }else if((cnt == '-' & side == 'outside') || (cnt == '+' & side == 'inside')){
    # extension towards lower genomic position
    end = start-seq(nchar(ins_seq)*1,nchar(ins_seq)*window_extend,1)+1
    if(cnt == '-' & side == 'outside'){
      # outside sequence starts one bp off of the breakend
      # inside is 'inclusive' of the breakend bp
      start = start - 1
      end = end - 1
    }
    align_scores = unlist(lapply(end,function(x) {
      refseq = getSeq(Hsapiens,chr,x,start)
      pairwiseAlignment(ins_seq,refseq,type = 'global-local',
                        substitutionMatrix = sub_mat,gapOpening = gap_open_pen,gapExtension = gap_ext_pen,scoreOnly = TRUE)
    }))
  }
  return(align_scores)
}

# given a SV breakend (location, +/-, window, alignment penalty parameters) and outside/inside, align insertions to k-2k bp of reference sequence and get (k+1) alignment scores
find_best_alignment_substring = function(ins_seq,window_extend,chr,cnt,side,ref_seq,gap_open_pen,gap_ext_pen,sub_mat){
  if((cnt == '+' & side == 'outside') || (cnt == '-' & side == 'inside')){
    # extension towards higher genomic position
    end = seq(nchar(ins_seq)*1,nchar(ins_seq)*window_extend,1)
    align_scores = unlist(lapply(substring(ref_seq,1,end),function(x) {
      pairwiseAlignment(ins_seq,DNAString(x),type = 'global-local',
                        substitutionMatrix = sub_mat,gapOpening = gap_open_pen,gapExtension = gap_ext_pen,scoreOnly = TRUE)
    }))
  }else if((cnt == '-' & side == 'outside') || (cnt == '+' & side == 'inside')){
    # extension towards lower genomic position
    end = seq(nchar(ins_seq)*(window_extend-1)+1,1,-1)
    align_scores = unlist(lapply(substring(ref_seq,end,nchar(ins_seq)*window_extend),function(x) {
      pairwiseAlignment(ins_seq,DNAString(x),type = 'global-local',
                        substitutionMatrix = sub_mat,gapOpening = gap_open_pen,gapExtension = gap_ext_pen,scoreOnly = TRUE)
    }))
  }
  return(align_scores)
}

# getting significance value from a alignment score matrix
alignment_sig = function(alignment_df,kmer){
  # kmer = 'TTTTGTTTAATTTAAATTAAAC'
  # alignment_df = fread('/Users/youyunzheng/Documents/HMS/PhD/beroukhimlab/broad_mount//youyun/nti/analysis_files/insertions/ins_align_total_02102315/TTTTGTTTAATTTAAATTAAAC/TTTTGTTTAATTTAAATTAAAC_alignment_quantile_out_og.tsv')
  # get the best alignment score (by quantile) for every breakend to get the adjusted p value
  best_alignment =
    data.table(melt(alignment_df,id.vars = 'SV_ID',variable.name = 'j',value.name = 'quantile'))[,k_mer := kmer][
    # setting j to the correct value (column 1 is supposed to be k by k, column 2 is supposed to be k by k+1. etc)
    ,.(j = as.numeric(gsub('V','',j)) + nchar(kmer) - 1,quantile,k_mer,
       # tagging max quantile value per breakend (SV_ID), max quantile in k by j at the breakend
       max_quantile = ifelse(quantile == max(quantile),TRUE,FALSE)),by = 'SV_ID'
  ][(max_quantile),
    # getting minimum avg and max offset with the max quantile value
    .(j_min = min(j),j_avg = mean(j),j_max = max(j))
    ,by = c('SV_ID','k_mer','max_quantile','quantile')]
  # best quantile per SV
  print(summary(best_alignment$quantile))
  # adjusting max qunatile
  best_alignment$adjusted_quantile = rank(best_alignment$quantile)/length(best_alignment$quantile)
  # adding significance info back into the total df
  best_alignment$significance = ifelse(best_alignment$adjusted_quantile > 0.95,TRUE,FALSE)
  return(best_alignment)
}

# align nearby function with multicore process
# given a SV call, generate alignment matricies for the SV call and their respective background
align_nearby_mc = function(ins_seq,window,insertion.sv.calls,intermediate_dir,gap_open,gap_epen,mismatch_pen,match_pen){
  window = as.numeric(window)
  mat =  nucleotideSubstitutionMatrix(match = match_pen, mismatch = mismatch_pen, type = 'DNA')[c('A','G','C','T','N'),c('A','G','C','T','N')]
  
  # creating intermeadiate directories to store all alignment matricies in
  intermediate_ins_dir = paste0(intermediate_dir,'/',ins_seq)
  dir.create(intermediate_ins_dir,showWarnings = TRUE)
  print(paste0('Creating ins directory for ',ins_seq,': ',intermediate_ins_dir))
  ins_alignment_scores = mclapply(X = c(1:nrow(insertion.sv.calls)),FUN = function(x) {
    each_call = as.character(unlist(insertion.sv.calls[x,]))
    chr = paste0('chr',each_call[which(colnames(insertion.sv.calls) == 'seqnames')])
    start = as.numeric(each_call[which(colnames(insertion.sv.calls) == 'start')])
    cnt = each_call[which(colnames(insertion.sv.calls) == 'cnt_type')]
    out_refseq = each_call[which(colnames(insertion.sv.calls) == 'outside_ref')]
    in_refseq = each_call[which(colnames(insertion.sv.calls) == 'inside_ref')]
    # return(list(find_best_alignment(DNAString(ins_seq),window,chr,start,cnt,'outside',gap_open,gap_epen,mat),
    #             find_best_alignment(DNAString(ins_seq),window,chr,start,cnt,'inside',gap_open,gap_epen,mat),
    #             find_best_alignment(reverseComplement(DNAString(ins_seq)),window,chr,start,cnt,'outside',gap_open,gap_epen,mat),
    #             find_best_alignment(reverseComplement(DNAString(ins_seq)),window,chr,start,cnt,'inside',gap_open,gap_epen,mat)))
    return(list(find_best_alignment_substring(DNAString(ins_seq),window,chr,cnt,'outside',out_refseq,gap_open,gap_epen,mat),
                find_best_alignment_substring(DNAString(ins_seq),window,chr,cnt,'inside',in_refseq,gap_open,gap_epen,mat),
                find_best_alignment_substring(reverseComplement(DNAString(ins_seq)),window,chr,cnt,'outside',out_refseq,gap_open,gap_epen,mat),
                find_best_alignment_substring(reverseComplement(DNAString(ins_seq)),window,chr,cnt,'inside',in_refseq,gap_open,gap_epen,mat)))
  },mc.cores = 4,mc.preschedule = TRUE,mc.set.seed = 55555)
  
  # Insertion alignment score matrix ----------
  ins_alignment_scores_out_og = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[1]]))
  ins_alignment_scores_in_og = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[2]]))
  ins_alignment_scores_out_rc = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[3]]))
  ins_alignment_scores_in_rc = do.call('rbind',lapply(ins_alignment_scores,function(x) x[[4]]))

  write.table(ins_alignment_scores_out_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_out_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_rc.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_scores_in_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_rc.tsv'),sep = '\t',row.names = FALSE)
  
  # column wise quantile matrix ----------
  # ins_alignment_scores_out_og <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_og.tsv'),sep = '\t')
  # ins_alignment_scores_in_og <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_og.tsv'),sep = '\t')
  # ins_alignment_scores_out_rc <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_out_rc.tsv'),sep = '\t')
  # ins_alignment_scores_in_rc <- fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_score_in_rc.tsv'),sep = '\t')
  
  ins_alignment_quantile_out_og = data.table(apply(ins_alignment_scores_out_og,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_in_og = data.table(apply(ins_alignment_scores_in_og,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_out_rc = data.table(apply(ins_alignment_scores_out_rc,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_in_rc = data.table(apply(ins_alignment_scores_in_rc,2,function(x) rank(x)/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  
  ins_alignment_quantile_out_og_max = data.table(apply(ins_alignment_scores_out_og,2,function(x) rank(x,ties.method = 'max')/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_in_og_max = data.table(apply(ins_alignment_scores_in_og,2,function(x) rank(x,ties.method = 'max')/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_out_rc_max = data.table(apply(ins_alignment_scores_out_rc,2,function(x) rank(x,ties.method = 'max')/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]
  ins_alignment_quantile_in_rc_max = data.table(apply(ins_alignment_scores_in_rc,2,function(x) rank(x,ties.method = 'max')/length(x)))[,SV_ID:= insertion.sv.calls$SV_ID]

  write.table(ins_alignment_quantile_out_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_in_og,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_out_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_rc.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_in_rc,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_rc.tsv'),sep = '\t',row.names = FALSE)

  write.table(ins_alignment_quantile_out_og_max,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_maxtie_out_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_in_og_max,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_maxtie_in_og.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_out_rc_max,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_maxtie_out_rc.tsv'),sep = '\t',row.names = FALSE)
  write.table(ins_alignment_quantile_in_rc_max,paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_maxtie_in_rc.tsv'),sep = '\t',row.names = FALSE)
  
  # row max quantile value calculation ----------
  # ins_alignment_quantile_out_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_og.tsv'),sep = '\t')
  # ins_alignment_quantile_in_og = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_og.tsv'),sep = '\t')
  # ins_alignment_quantile_out_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_out_rc.tsv'),sep = '\t')
  # ins_alignment_quantile_in_rc = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_alignment_quantile_in_rc.tsv'),sep = '\t')

  out.ins.match = alignment_sig(ins_alignment_quantile_out_og,ins_seq)
  in.ins.match = alignment_sig(ins_alignment_quantile_in_og,ins_seq)
  out.ins.rc.match = alignment_sig(ins_alignment_quantile_out_rc,ins_seq)
  in.ins.rc.match = alignment_sig(ins_alignment_quantile_in_rc,ins_seq)

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
  window = as.numeric(opt$window)
  insertion.sv.calls = fread(opt$data)
  intermediate_dir = opt$output
  
  if(!is.numeric(window) || is.null(ins) || nrow(insertion.sv.calls) == 0 || intermediate_dir == ''){
    print('check inputs!!')
  }
  
  # insertion.sv.calls.subset = insertion.sv.calls[c(1:1000)]
  # insertion.sv.calls[,c('outside_ref','inside_ref') := list(
  #   find_surrounding_seq(nchar(ins)*window,paste0('chr',seqnames),start,cnt_type,'outside'),
  #   find_surrounding_seq(nchar(ins)*window,paste0('chr',seqnames),start,cnt_type,'inside')
  # ),by = .(SV_ID)]
    
  print(paste0('Aligning the insertion sequence [',ins,'] with a ',window,'X window using all breakends in file [',opt$data,'] and writing to [',intermediate_dir,']'))
  system.time(unlist(align_nearby_mc(ins,window,insertion.sv.calls,intermediate_dir,7,1,1,3)))
  
}