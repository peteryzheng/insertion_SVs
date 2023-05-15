library(VennDiagram)
library(data.table)
library(UpSetR)

# local vs UGER
if(Sys.getenv("LOGNAME") == 'youyunzheng'){
  workdir = '~/Documents/HMS/PhD/beroukhimlab/broad_mount/'
}else{
  workdir = '/xchip/beroukhimlab/'
}

intermediate_dir = paste0(workdir,'/youyun/nti/analysis_files/insertions')
dataset = 'DIPG'
dataset = 'TCGA'
if(dataset == 'DIPG'){
  # DIPG
  insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_050916.tsv'))
  intermediate_dir = paste0(intermediate_dir,'/ins_align_total_05102316')
}else{
  # TCGA
  # insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_050918.tsv'))
  # intermediate_dir = paste0(intermediate_dir,'/ins_align_total_0510231628')
  insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_030900.tsv'))
  intermediate_dir = paste0(intermediate_dir,'/ins_align_total_03092312')
  
}

insertion.sv.calls.subset = insertion.sv.calls[ins_len <= 30 & ins_len >= 6]

# single breakend matching calling ---------- 
match.results = data.table(t(apply(insertion.sv.calls.subset,1,function(x){
  # x = as.character(insertion.sv.calls.subset[2,])
  ins_seq = x[which(colnames(insertion.sv.calls) == 'ins_seq')]
  current_SV_ID = x[which(colnames(insertion.sv.calls) == 'SV_ID')]
  intermediate_ins_dir = paste0(intermediate_dir,'/',ins_seq)

  out_og_breakend_results = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_out_og_sig_breakends.tsv'))
  in_og_breakend_results = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_in_og_sig_breakends.tsv'))
  out_rc_breakend_results = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_out_rc_sig_breakends.tsv'))
  in_rc_breakend_results = fread(paste0(intermediate_ins_dir,'/',ins_seq,'_in_rc_sig_breakends.tsv'))

  # Checking if out SV_ID is in the significant match breakend list in the output in the 4 possible configurations
  out.ins.match = c(out_og_breakend_results[SV_ID == current_SV_ID]$quantile,
                    out_og_breakend_results[SV_ID == current_SV_ID]$j_min,
                    ifelse(current_SV_ID %in% out_og_breakend_results[(significance)]$SV_ID,1,0))
  in.ins.match = c(in_og_breakend_results[SV_ID == current_SV_ID]$quantile,
                   in_og_breakend_results[SV_ID == current_SV_ID]$j_min,
                   ifelse(current_SV_ID %in% in_og_breakend_results[(significance)]$SV_ID,1,0))
  out.ins.rc.match = c(out_rc_breakend_results[SV_ID == current_SV_ID]$quantile,
                       out_rc_breakend_results[SV_ID == current_SV_ID]$j_min,
                       ifelse(current_SV_ID %in% out_rc_breakend_results[(significance)]$SV_ID,1,0))
  in.ins.rc.match = c(in_rc_breakend_results[SV_ID == current_SV_ID]$quantile,
                      in_rc_breakend_results[SV_ID == current_SV_ID]$j_min,
                      ifelse(current_SV_ID %in% in_rc_breakend_results[(significance)]$SV_ID,1,0))

  return(c(current_SV_ID,out.ins.match,in.ins.match,out.ins.rc.match,in.ins.rc.match))
})))
# merge the result back into the original sv file
colnames(match.results) = c('SV_ID','outside_ins_match_quantile','outside_ins_match_j_max_quantile','outside_ins_match',
                            'inside_ins_match_quantile','inside_ins_match_j_max_quantile','inside_ins_match',
                            'outside_ins_rc_match_quantile','outside_ins_rc_match_j_max_quantile','outside_ins_rc_match',
                            'inside_ins_rc_match_quantile','inside_ins_rc_match_j_max_quantile','inside_ins_rc_match')
match.results = match.results[,lapply(.SD,as.numeric),by=SV_ID]
insertion.sv.calls.aligned = data.table(merge(insertion.sv.calls.subset,match.results,by = 'SV_ID'))

insertion.sv.calls.aligned$seqnames = factor(insertion.sv.calls.aligned$seqnames,levels = c(seq(1,22,1),'X'))
insertion.sv.calls.aligned[,c('chr_order','SV_config_combo') := list(ifelse(seqnames[1] <= seqnames[2],TRUE,FALSE),paste0(cnt_type,collapse = '')),
                           by = .(SV_ID = gsub(':.*','',ID),Sample)]
write.table(insertion.sv.calls.aligned,paste0(workdir,'youyun/nti/analysis_files/',dataset,'_insertions_w_sig_alignment_',format(Sys.time(), "%m%d%H"),'.tsv'),
            sep = '\t',row.names = FALSE)
# making a upset plot on the breakend matching result
pdf(paste0(workdir,'/youyun/nti/analysis_files/',dataset,'_alignment_upset_plot_',format(Sys.time(), "%m%d%y%H"),'.pdf'),width=8,height=8)
upset(insertion.sv.calls.aligned,sets = c('outside_ins_match','inside_ins_match','outside_ins_rc_match','inside_ins_rc_match'),
      sets.bar.color = "#56B4E9")
dev.off()
# insertion.sv.calls.aligned = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_w_sig_alignment_030916.tsv'))
# write.table(insertion.sv.calls.aligned,paste0(workdir,'youyun/nti/analysis_files/insertions_w_sig_alignment_',format(Sys.time(), "%m%d%H"),'.csv'),
#             sep = ',',row.names = FALSE)

# function to find convolution p value cut off ----------
convolve_sig = function(ins_seq,intermediate_dir,side1,side2){
  # ins_seq = 'ACACCTTGC'
  # side1 = 'out_og'
  # side2 = 'out_og'
  # we are calculating convolution background by convolving the two sides (i.e. outog with outog) within a SV breakend pair
  side1_quants = fread(paste0(intermediate_dir,'/',ins_seq,'/',ins_seq,'_',side1,'_sig_breakends.tsv'))[,.(SV_ID,side1_quantile = quantile)]
  side2_quants = fread(paste0(intermediate_dir,'/',ins_seq,'/',ins_seq,'_',side2,'_sig_breakends.tsv'))[,.(SV_ID,side2_quantile = quantile)]
  merge_quants = merge(side1_quants,side2_quants,by = 'SV_ID')
  # averaging between side1 quantile of breakend 1 and side 2 quantile of breakend 2, then side 1 quantile of breakend 2 and side 2 quantile of breakend 1
  conv_avg_quants = merge_quants[,.(mean1 = (side1_quantile[1]+side2_quantile[2])/2,mean2 = (side1_quantile[2]+side2_quantile[1])/2),
                                 # per breakend pair of SV
                                 by = gsub('__.*.__|:1|:2','',SV_ID)]
  # min value between side1 quantile of breakend 1 and side 2 quantile of breakend 2, then side 1 quantile of breakend 2 and side 2 quantile of breakend 1
  conv_min_quants = merge_quants[,.(min1 = min(side1_quantile[1],side2_quantile[2]),min2 = min(side1_quantile[2],side2_quantile[1])),
                                 # per breakend pair of SV
                                 by = gsub('__.*.__|:1|:2','',SV_ID)]
  c(cutoff_avg = as.numeric(quantile(c(conv_avg_quants$mean1,conv_avg_quants$mean2),0.95)),
    cutoff_min = as.numeric(quantile(c(conv_min_quants$min1,conv_min_quants$min2),0.95)))
}

# outside + outside -------------------------------------------------------
# # convolution implementation testing
# library(ggplot2)
# ins_seq = 'ACACCTTGC'
# side1 = 'out_og'
# side2 = 'out_rc'
# side1_probs = fread(paste0(intermediate_dir,'/',ins_seq,'/',ins_seq,'_',side1,'_sig_breakends.tsv'))$quantile
# side2_probs = fread(paste0(intermediate_dir,'/',ins_seq,'/',ins_seq,'_',side2,'_sig_breakends.tsv'))$quantile
# conv_sum_probs = side1_probs + side2_probs
# conv_min_probs = mapply(min,side1_probs,side2_probs)
# summary(side1_probs)
# summary(side2_probs)
# summary(conv_sum_probs)
# summary(conv_min_probs)
# ggplot(melt(data.table(SV = c(1:length(side1_probs)),side1_probs,side2_probs,conv_sum_probs/2,conv_min_probs),id.vars = 'SV')) + 
#   geom_density(aes(x = value,fill = variable),alpha = 0.2)
# 
# outside outside w convolution
# first for each SV pair, get the background convolution cutoff value for each respective SV type
insertion.sv.calls.aligned[,c('cutoff_avg','cutoff_min') := as.list(
  ifelse(SV_config_combo %in% c('-+','+-'),
         # for -/+ or +/- we use outside og outside og as background for convolution
         convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_og'),
         # for -/- or +/+ we use outside og outside rc as background for convolution
         convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_rc'))
),by = .(SV_ID = gsub(':.*','',ID),Sample)]

insertion.sv.calls.aligned[
  ,c('outout_ins_avg','outout_ins_min') := list(
    ifelse(SV_config_combo %in% c('-+','+-') &
             # if +/- or -/+, then outside og and outside og convolution
             (mean(outside_ins_match_quantile) > unique(cutoff_avg)),
           1,ifelse(SV_config_combo %in% c('--','++') &
                         # if +/+ or -/-, then outside og and outside rc convolution
                         ((mean(c(outside_ins_match_quantile[1],outside_ins_rc_match_quantile[2])) > unique(cutoff_avg)) |
                            (mean(c(outside_ins_match_quantile[2],outside_ins_rc_match_quantile[1])) > unique(cutoff_avg))),
                       1,0)),
    ifelse(SV_config_combo %in% c('-+','+-') &
             # if +/- or -/+, then outside og and outside og convolution
             (min(outside_ins_match_quantile) > unique(cutoff_min)),
           1,ifelse(SV_config_combo %in% c('--','++') &
                         # if +/+ or -/-, then outside og and outside rc convolution
                         ((min(c(outside_ins_match_quantile[1],outside_ins_rc_match_quantile[2])) > unique(cutoff_min)) |
                            (min(c(outside_ins_match_quantile[2],outside_ins_rc_match_quantile[1])) > unique(cutoff_min))),
                       1,0))
  ),by = .(SV_ID = gsub(':.*','',ID),Sample)
]
table(insertion.sv.calls.aligned$outout_ins_avg)
table(insertion.sv.calls.aligned$outout_ins_min)
write.table(insertion.sv.calls.aligned[outout_ins_avg | outout_ins_min],
            paste0(workdir,'youyun/nti/analysis_files/',dataset,'_outside_outside_aligned_',format(Sys.time(), "%m%d%H"),'.tsv'),
            sep = '\t',row.names = FALSE)

upset_df = dcast(melt(insertion.sv.calls.aligned[,.(SV_ID = paste0(Sample,'__',gsub(':.*','',ID)),SV_config_combo,outout_ins_avg,outout_ins_min,
                                                    og_rc_match = paste0(gsub('.*:','',ID),':',outside_ins_match,'_',outside_ins_rc_match))],
                      id.vars = 'SV_ID')[,value := paste0(variable,': ',value)],SV_ID~value,fill = 0,fun = function(x){ifelse(length(x) > 0,1,0)})
pdf(paste0(workdir,'/youyun/nti/analysis_files/',dataset,'_alignment_outside_outside_upset_plot_',format(Sys.time(), "%m%d%y%H"),'.pdf'),width=8,height=8)
upset(upset_df,sets = c('SV_config_combo: ++','SV_config_combo: --','SV_config_combo: +-','SV_config_combo: -+'),
      sets.bar.color = "#56B4E9")
upset(upset_df,sets = c('SV_config_combo: ++','SV_config_combo: --','SV_config_combo: +-','SV_config_combo: -+',
                        "outout_ins_avg: 0","outout_ins_avg: 1","outout_ins_min: 0","outout_ins_min: 1"),
      sets.bar.color = "#56B4E9",keep.order = TRUE)
upset(upset_df,sets = c('SV_config_combo: ++','SV_config_combo: --','SV_config_combo: +-','SV_config_combo: -+',
                        "og_rc_match: 1:0_0","og_rc_match: 1:0_1","og_rc_match: 1:1_0","og_rc_match: 1:1_1",
                        "og_rc_match: 2:0_0","og_rc_match: 2:0_1","og_rc_match: 2:1_0","og_rc_match: 2:1_1"),
      sets.bar.color = "#56B4E9",keep.order = TRUE,order.by = 'degree',nintersects = NA,decreasing = TRUE)
upset(upset_df,sets = c('SV_config_combo: ++','SV_config_combo: --','SV_config_combo: +-','SV_config_combo: -+',
                        "outout_ins_avg: 0","outout_ins_avg: 1","outout_ins_min: 0","outout_ins_min: 1",
                        "og_rc_match: 1:0_0","og_rc_match: 1:0_1","og_rc_match: 1:1_0","og_rc_match: 1:1_1",
                        "og_rc_match: 2:0_0","og_rc_match: 2:0_1","og_rc_match: 2:1_0","og_rc_match: 2:1_1"),
      sets.bar.color = "#56B4E9",keep.order = TRUE,order.by = 'degree',nintersects = NA,decreasing = TRUE)
dev.off()

# old outside outside calling 
# insertion.sv.calls.aligned[,outout_ins := ifelse(
#   # if connection is +- or -+ we want to see either (both outside ins match)
#   paste0(cnt_type,collapse = '') %in% c('+-','-+') &  
#     (all(outside_ins_match>0) | all(outside_ins_rc_match>0)),
#   TRUE,ifelse(
#     # if connection is ++ or -- we want to see either (outside ins + outside rc ins) or (outside rc ins + outside ins)
#     paste0(cnt_type,collapse = '') %in% c('++','--') &  
#       (all(c(outside_ins_match[1],outside_ins_rc_match[2])>0) | 
#          all(c(outside_ins_match[2],outside_ins_rc_match[1])>0)),
#     TRUE,FALSE
#   )),
#   by = .(SV_ID = gsub(':.*','',ID),Sample)] 
# table(insertion.sv.calls.aligned$outout_ins)
# upset(dcast(insertion.sv.calls.aligned[,.(SV_ID,SV_config_combo,outout_ins = ifelse(outout_ins,1,0))],
#             formula = SV_ID+outout_ins~SV_config_combo,fill = 0),
#       sets = c('++','--','+-','-+','outout_ins'), 
#       sets.bar.color = "#56B4E9",empty.intersections = TRUE)

# # getting all the correct outside inside matches ----------
# insertion.sv.calls.aligned[
#   ,c('outout_ins_avg','outout_ins_min') := list(
#     ifelse(SV_config_combo %in% c('-+','+-') & 
#              # if +/- or -/+, then outside og and outside og convolution
#              (mean(outside_ins_match_quantile) > 
#                 convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_og')["cutoff_avg"]),
#            TRUE,ifelse(SV_config_combo %in% c('--','++') & 
#                          # if +/+ or -/-, then outside og and outside rc convolution
#                          ((mean(c(outside_ins_match_quantile[1],outside_ins_rc_match_quantile[2])) > 
#                              convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_rc')["cutoff_avg"]) | 
#                             (mean(c(outside_ins_match_quantile[2],outside_ins_rc_match_quantile[1])) > 
#                                convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_rc')["cutoff_avg"])),
#                        TRUE,FALSE)),
#     ifelse(SV_config_combo %in% c('-+','+-') & 
#              # if +/- or -/+, then outside og and outside og convolution
#              (min(outside_ins_match_quantile) > 
#                 convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_og')["cutoff_min"]),
#            TRUE,ifelse(SV_config_combo %in% c('--','++') & 
#                          # if +/+ or -/-, then outside og and outside rc convolution
#                          ((min(c(outside_ins_match_quantile[1],outside_ins_rc_match_quantile[2])) > 
#                              convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_rc')["cutoff_min"]) | 
#                             (min(c(outside_ins_match_quantile[2],outside_ins_rc_match_quantile[1])) > 
#                                convolve_sig(unique(ins_seq),intermediate_dir,'out_og','out_rc')["cutoff_min"])),
#                        TRUE,FALSE))
#   ),by = .(SV_ID = gsub(':.*','',ID),Sample)
# ]
# 
# old outside inside calling
# insertion.sv.calls.aligned[,outin_ins := ifelse(
#   # if connection is +- or -+ we want to see either (outside ins + inside ins match) or (inside ins + outside ins match)
#   paste0(cnt_type,collapse = '') %in% c('+-','-+') &
#     (all(c(outside_ins_match[1],inside_ins_match[2])>0) |
#        all(c(outside_ins_match[2],inside_ins_match[1])>0)),
#   TRUE,ifelse(
#     # if connection is ++ or -- we want to see either (outside ins + inside rc ins) or (outside rc ins + inside ins)
#     # or (inside ins + outside rc ins) or (inside rc ins + outside ins)
#     paste0(cnt_type,collapse = '') %in% c('++','--') &
#       (all(c(outside_ins_match[1],inside_ins_rc_match[2])>0) |
#          all(c(outside_ins_rc_match[1],inside_ins_match[2])>0) |
#          all(c(outside_ins_match[2],inside_ins_rc_match[1])>0) |
#          all(c(outside_ins_rc_match[2],inside_ins_match[1])>0)
#       ),
#     TRUE,FALSE
#   )),
#   by = .(SV_ID = gsub(':.*','',ID),Sample)]
# table(insertion.sv.calls.aligned$outin_ins)
# write.table(insertion.sv.calls.aligned[outin_ins & ins_len > 4],
#             paste0(workdir,'youyun/nti/analysis_files/DIPG_outside_inside_aligned_',format(Sys.time(), "%m%d%H"),'.tsv'),
#             sep = '\t',row.names = FALSE)
# 
# # getting all the correct inside inside matches ----------
# insertion.sv.calls.aligned[,inin_ins := ifelse(
#   # if connection is +- or -+ we want to see either (both inside ins match) or (both outside rc ins match)
#   paste0(cnt_type,collapse = '') %in% c('+-','-+') &  
#     (all(inside_ins_match>0)| all(inside_ins_rc_match > 0)),
#   TRUE,ifelse(
#     # if connection is ++ or -- we want to see either (outside ins + outside rc ins) or (outside rc ins + outside ins)
#     paste0(cnt_type,collapse = '') %in% c('++','--') &  
#       (all(c(inside_ins_match[1],inside_ins_rc_match[2])>0) | 
#          all(c(inside_ins_match[2],inside_ins_rc_match[1])>0)),
#     TRUE,FALSE
#   )),
#   by = .(SV_ID = gsub(':.*','',ID),Sample)] 
# table(insertion.sv.calls.aligned$inin_ins)
# write.table(insertion.sv.calls.aligned[inin_ins & ins_len > 4],
#             paste0(workdir,'youyun/nti/analysis_files/DIPG_inside_inside_aligned_',format(Sys.time(), "%m%d%H"),'.tsv'),
#             sep = '\t',row.names = FALSE)
# 
# venn.diagram(
#   x = list(insertion.sv.calls.aligned[outout_ins == TRUE]$SV_ID,insertion.sv.calls.aligned[outin_ins == TRUE]$SV_ID,insertion.sv.calls.aligned[inin_ins == TRUE]$SV_ID),
#   category.names = c('Outside outside','Outside inside','Inside inside'),
#   filename = paste0(workdir,'youyun/nti/analysis_files/DIPG_configuration_',format(Sys.time(), "%m%d%H"),'.png'),output = TRUE,
#   height = 480,width = 480,resolution = 300,compression = "lzw",col=c("#440154ff", '#21908dff', '#fde725ff'),
#   # Circles
#   lwd = 2,lty = 'blank',fill = c("#440154ff", '#21908dff', '#fde725ff'),
#   # Numbers
#   cex = .6,fontface = "bold",fontfamily = "sans",
#   # set names
#   cat.cex = 0.3,cat.fontface = "bold",cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
#   cat.dist = c(0.055, 0.055, 0.085),cat.fontfamily = "sans",rotation = 1
# )
# 
