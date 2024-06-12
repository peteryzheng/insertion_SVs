FROM continuumio/anaconda3:2024.02-1@sha256:a433c3912e65f020377c2b553c26860fd7e5c44ff854273cc37a43ad6dc48e99

RUN conda config --add channels defaults  && \
    conda config --add channels bioconda  && \
    conda config --add channels conda-forge  && \
    conda config --set channel_priority strict  && \
    conda create -n insertions -y \
    r-base r-data.table r-optparse \
    r-ggplot2 r-rcolorbrewer r-ggrepel \
    r-upsetr r-parallelly r-stringr \
    git \
    bioconductor-bsgenome \
    bioconductor-bsgenome.hsapiens.ucsc.hg19 \
    bioconductor-bsgenome.hsapiens.ucsc.hg38 && \
    conda run -n insertions git clone https://github.com/peteryzheng/insertion_SVs.git

