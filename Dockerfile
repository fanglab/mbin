# Install base with hdf5, python, R to ease installation of Pacbio tools and libraries
FROM mpkocher/docker-pacbiobase

################## METADATA ######################
LABEL maintainer="john.beaulaurier@gmail.com"
LABEL base.image="mpkocher/docker-pacbiobase"
LABEL version="1"
LABEL software="mbin"
LABEL software.version="1.1.1"
LABEL description="Tool for methylation profiling and binning of metagenomic sequences"
LABEL website="https://github.com/fanglab/mbin"
LABEL documentation="https://mbin.readthedocs.io/en/latest/"
LABEL license="https://github.com/fanglab/mbin/blob/master/LICENSE"
LABEL tags="Genomics,Methylation"

################# INSTALLATION #####################

# Install pbcore
RUN python -m pip install --upgrade pip setuptools wheel
WORKDIR /tmp
RUN git clone https://github.com/PacificBiosciences/pbcore.git
WORKDIR /tmp/pbcore
RUN pip install .

# Install scipy, matplotlib, pysam, and biopython
RUN pip install scipy matplotlib pysam biopython

# Install mbin
RUN pip install mbin

# Install bh-tsne
WORKDIR /tmp
RUN git clone https://github.com/lvdmaaten/bhtsne.git
WORKDIR /tmp/bhtsne 
RUN g++ sptree.cpp tsne.cpp tsne_main.cpp -o bh_tsne -O2 
ENV PATH=$PATH:/tmp/bhtsne

# Set mbin working directory
WORKDIR /mbin
