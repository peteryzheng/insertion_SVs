library(VennDiagram)
library(data.table)

# local vs UGER
if(Sys.getenv("LOGNAME") == 'youyunzheng'){
  workdir = '~/Documents/HMS/PhD/beroukhimlab/broad_mount/'
}else{
  workdir = '/xchip/beroukhimlab/'
}

intermediate_dir = paste0(workdir,'/youyun/nti/analysis_files/insertions')
# create intermediate directory to store all intermediate alignment results for RAM efficiency
intermediate_dir = paste0(intermediate_dir,'/ins_align_total_03012300')
insertion.sv.calls = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_020716.tsv'))

match.results = t(apply(insertion.sv.calls[ins_len <= 30 & ins_len >= 6],1,function(x){
  # x = as.character(insertion.sv.calls[ins_len <= 30 & ins_len >= 6][1,])
  ins_seq = x[which(colnames(insertion.sv.calls) == 'ins_seq')]
  SV_ID = x[which(colnames(insertion.sv.calls) == 'SV_ID')]
  intermediate_ins_dir = paste0(intermediate_dir,'/',ins_seq)
  
  # Checking if out SV_ID is in the significant match breakend list in the output in the 4 possible configurations
  out.ins.match = ifelse(SV_ID %in% fread(paste0(intermediate_ins_dir,'/',ins_seq,'_out_og_sig_breakends.tsv'))[(significance)]$SV_ID,1,0)
  in.ins.match = ifelse(SV_ID %in% fread(paste0(intermediate_ins_dir,'/',ins_seq,'_in_og_sig_breakends.tsv'))[(significance)]$SV_ID,1,0)
  out.ins.rc.match = ifelse(SV_ID %in% fread(paste0(intermediate_ins_dir,'/',ins_seq,'_out_rc_sig_breakends.tsv'))[(significance)]$SV_ID,1,0)
  in.ins.rc.match = ifelse(SV_ID %in% fread(paste0(intermediate_ins_dir,'/',ins_seq,'_in_rc_sig_breakends.tsv'))[(significance)]$SV_ID,1,0)
  
  return(c(SV_ID,out.ins.match,in.ins.match,out.ins.rc.match,in.ins.rc.match))
}))

colnames(match.results) = c('SV_ID','outside_ins_match','inside_ins_match','outside_ins_rc_match','inside_ins_rc_match')
insertion.sv.calls.aligned = data.table(merge(insertion.sv.calls[ins_len <= 30 & ins_len >= 6],match.results,by = 'SV_ID'))

write.table(insertion.sv.calls.aligned,paste0(workdir,'youyun/nti/analysis_files/insertions_w_sig_alignment_',format(Sys.time(), "%m%d%H"),'.tsv'),
            sep = '\t',row.names = FALSE)

pdf(paste0(workdir,'/youyun/nti/analysis_files/DIPG_alignment_upset_plot_',format(Sys.time(), "%m%d%y%H"),'.pdf'),width=8,height=8)
upset(insertion.sv.calls.aligned,sets = c('outside_ins_match','inside_ins_match','outside_ins_rc_match','inside_ins_rc_match'), 
      sets.bar.color = "#56B4E9",empty.intersections = 'on')
dev.off()

# insertion.sv.calls.aligned = fread(paste0(workdir,'youyun/nti/analysis_files/insertions_w_sig_alignment_022811.tsv'))
# write.table(insertion.sv.calls.aligned,paste0(workdir,'youyun/nti/analysis_files/insertions_w_sig_alignment_',format(Sys.time(), "%m%d%H"),'.csv'),
#             sep = ',',row.names = FALSE)


# outside + outside -------------------------------------------------------

# focusing on insertions for now
any(is.na(insertion.sv.calls.aligned[,c(20:23)]))
colnames(insertion.sv.calls.aligned)[c(20:23)]

# getting all the correct outside outside matches
insertion.sv.calls.aligned[,outout_ins := ifelse(
  # if connection is +- or -+ we want to see either (both outside ins match)
  paste0(cnt_type,collapse = '') %in% c('+-','-+') &  
    (all(outside_ins_match>0) | all(outside_ins_rc_match>0)),
  TRUE,ifelse(
    # if connection is ++ or -- we want to see either (outside ins + outside rc ins) or (outside rc ins + outside ins)
    paste0(cnt_type,collapse = '') %in% c('++','--') &  
      (all(c(outside_ins_match[1],outside_ins_rc_match[2])>0) | 
         all(c(outside_ins_match[2],outside_ins_rc_match[1])>0)),
    TRUE,FALSE
  )),
  by = .(SV_ID = gsub(':.*','',ID),Sample)] 
table(insertion.sv.calls.aligned$outout_ins)

# getting all the correct inside inside matches
insertion.sv.calls.aligned[,inin_ins := ifelse(
  # if connection is +- or -+ we want to see either (both inside ins match) or (both outside rc ins match)
  paste0(cnt_type,collapse = '') %in% c('+-','-+') &  
    (all(inside_ins_match>0)| all(inside_ins_rc_match > 0)),
  TRUE,ifelse(
    # if connection is ++ or -- we want to see either (outside ins + outside rc ins) or (outside rc ins + outside ins)
    paste0(cnt_type,collapse = '') %in% c('++','--') &  
      (all(c(inside_ins_match[1],inside_ins_rc_match[2])>0) | 
         all(c(inside_ins_match[2],inside_ins_rc_match[1])>0)),
    TRUE,FALSE
  )),
  by = .(SV_ID = gsub(':.*','',ID),Sample)] 
table(insertion.sv.calls.aligned$inin_ins)

# getting all the correct outside inside matches
insertion.sv.calls.aligned[,outin_ins := ifelse(
  # if connection is +- or -+ we want to see either (outside ins + inside ins match) or (inside ins + outside ins match)
  paste0(cnt_type,collapse = '') %in% c('+-','-+') &  
    (all(c(outside_ins_match[1],inside_ins_match[2])>0) | 
       all(c(outside_ins_match[2],inside_ins_match[1])>0)),
  TRUE,ifelse(
    # if connection is ++ or -- we want to see either (outside ins + inside rc ins) or (outside rc ins + inside ins) 
    # or (inside ins + outside rc ins) or (inside rc ins + outside ins)
    paste0(cnt_type,collapse = '') %in% c('++','--') &  
      (all(c(outside_ins_match[1],inside_ins_rc_match[2])>0) | 
         all(c(outside_ins_rc_match[1],inside_ins_match[2])>0) | 
         all(c(outside_ins_match[2],inside_ins_rc_match[1])>0) | 
         all(c(outside_ins_rc_match[2],inside_ins_match[1])>0)
      ),
    TRUE,FALSE
  )),
  by = .(SV_ID = gsub(':.*','',ID),Sample)] 
table(insertion.sv.calls.aligned$outin_ins)


write.table(insertion.sv.calls.aligned[outout_ins & ins_len > 4],
            paste0(workdir,'youyun/nti/analysis_files/DIPG_outside_outside_aligned_',format(Sys.time(), "%m%d%H"),'.tsv'),
            sep = '\t',row.names = FALSE)
write.table(insertion.sv.calls.aligned[inin_ins & ins_len > 4],
            paste0(workdir,'youyun/nti/analysis_files/DIPG_inside_inside_aligned_',format(Sys.time(), "%m%d%H"),'.tsv'),
            sep = '\t',row.names = FALSE)
write.table(insertion.sv.calls.aligned[outin_ins & ins_len > 4],
            paste0(workdir,'youyun/nti/analysis_files/DIPG_outside_inside_aligned_',format(Sys.time(), "%m%d%H"),'.tsv'),
            sep = '\t',row.names = FALSE)


venn.diagram(
  x = list(insertion.sv.calls.aligned[outout_ins == TRUE]$SV_ID,insertion.sv.calls.aligned[outin_ins == TRUE]$SV_ID,insertion.sv.calls.aligned[inin_ins == TRUE]$SV_ID),
  category.names = c('Outside outside','Outside inside','Inside inside'),
  filename = paste0(workdir,'youyun/nti/analysis_files/DIPG_configuration_',format(Sys.time(), "%m%d%H"),'.png'),output = TRUE,
  height = 480,width = 480,resolution = 300,compression = "lzw",col=c("#440154ff", '#21908dff', '#fde725ff'),
  # Circles
  lwd = 2,lty = 'blank',fill = c("#440154ff", '#21908dff', '#fde725ff'),
  # Numbers
  cex = .6,fontface = "bold",fontfamily = "sans",
  # set names
  cat.cex = 0.3,cat.fontface = "bold",cat.default.pos = "outer",cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),cat.fontfamily = "sans",rotation = 1
)

