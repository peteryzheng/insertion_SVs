library(data.table)
library(ggplot2)
library(RColorBrewer)
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

# DATA ----------
# PCAWG path
# sv_files = list.files(paste0(workdir,'/siyun/data/insertions/pcawg/'),
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
# filter for minimal qual score and samples with tumor in normal contamination
sv.calls = sv.calls[as.numeric(QUAL) > 15 & !Sample %in% c('DIPG19_TOR_pair', 'DIPG58_TOR_pair', 
                                                           'SJHGG008_A_STJ_WGS_pair', 'SJHGG041_D_STJ_WGS_pair'),]
# filter for only the major chromosomes
sv.calls = sv.calls[seqnames %in% c(1:22,'X','Y'),]
sv.calls[,SV_ID := paste(Sample,seqnames,ID,sep = '__')]
# making sure every breakend will have its pair (both breakends are on chr 1-22)
sv.calls[,both_end_pass_filter := ifelse(.N == 2,TRUE,FALSE),,by = .(SV_pair_ID = gsub(':.*','',ID),Sample)] 
sv.calls = sv.calls[both_end_pass_filter == TRUE]
# filter for sv breakpoints with insertion sequences only (if one breakend has insertion in info, the paired breakend will also have insertion in info)
insertion.sv.calls = sv.calls[grepl('INSERTION',INFO),]

print(paste0('Total number of SV breakends (qual score <= 15 discarded): ',nrow(sv.calls)))
# -4 bc 4 samples have tumor in normal contamination
print(paste0('Total number of samples before insertion filter: ',length(sv_files)-4))
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

insertion.sv.calls[,c('outside_ref','inside_ref') := list(
  find_surrounding_seq(30*5,paste0('chr',seqnames),start,cnt_type,'outside'),
  find_surrounding_seq(30*5,paste0('chr',seqnames),start,cnt_type,'inside')
),by = .(SV_ID)]

ggplot(melt(merge(sv.calls[,.(total_SV = .N/2),Sample],insertion.sv.calls[,.(ins_SV = .N/2),Sample]),id.vars = 'Sample'),
       aes(x = Sample,y = value,fill = variable)) + 
  geom_bar(position = "dodge", stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1 , size = 5)) +
  scale_fill_brewer(palette = "Set1")+
  labs(x = "Sample", y = "Count") +
  ggtitle("Histogram of total_SV and ins_SV by Sample")
ggsave(paste0(workdir,'youyun/nti/analysis_files/SV_burden_',format(Sys.time(), "%m%d%H"),'.pdf'), 
       plot = last_plot(), device = "pdf")

sv.file.path = paste0(workdir,'youyun/nti/analysis_files/total_SVs_processed_',format(Sys.time(), "%m%d%H"),'.tsv')
ins.file.path = paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_',format(Sys.time(), "%m%d%H"),'.tsv')
ins.file.filter.hypermut.path = paste0(workdir,'youyun/nti/analysis_files/insertions_SVs_processed_filter_hypermut_',
                                       format(Sys.time(), "%m%d%H"),'.tsv')
write.table(sv.calls,sv.file.path,sep = '\t',row.names = FALSE)
write.table(insertion.sv.calls,ins.file.path,sep = '\t',row.names = FALSE)
write.table(insertion.sv.calls[Sample %in% sv.calls[,.N,Sample][N<= 500,]$Sample],
            ins.file.filter.hypermut.path,sep = '\t',row.names = FALSE)
print(paste0('Total number of samples after filter: ',
             length(unique(insertion.sv.calls[Sample %in% sv.calls[,.N,Sample][N<= 500,]$Sample]$Sample))))
print(paste0('output file here: ',ins.file.path))