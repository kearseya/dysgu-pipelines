FROM ubuntu:jammy

ENV DEBIAN_FRONTEND=noninteractive

# 3.10 not in noble
# RUN apt-get update && apt-get install -y --no-install-recommends software-properties-common
# RUN add-apt-repository ppa:deadsnakes/ppa

RUN apt-get update && apt-get install -y --no-install-recommends build-essential r-base r-cran-randomforest python3.10 python3-pip python3-setuptools python3-dev

RUN apt-get install -y --no-install-recommends git vim make

RUN apt-get install -y --no-install-recommends cython3 python3-pip

RUN apt-get install -y --no-install-recommends libghc-zlib-dev zlib1g-dev zlib1g libbz2-dev liblzma-dev libcurl4-openssl-dev libcrypto++-dev libssl-dev libdeflate-dev automake autoconf cmake 

RUN apt-get install -y --no-install-recommends bedtools poppler-utils imagemagick ghostscript

WORKDIR /app

# htslib
RUN git clone https://github.com/samtools/htslib.git; cd htslib; git submodule update --init --recursive; autoreconf -i; ./configure; make; make install; cd /app

# bcftools
RUN git clone https://github.com/samtools/bcftools.git; cd bcftools; autoreconf -i; ./configure; make; make install; cd /app

# samtools
RUN git clone https://github.com/samtools/samtools.git; cd samtools; autoheader; autoconf -Wno-syntax; ./configure; make; make install; cd /app

# bwa
RUN cd /app; git clone https://github.com/lh3/bwa.git; cd bwa; make; make install; cp bwa /usr/local/bin; cd /app

# dodi
RUN git clone https://github.com/kcleal/dodi.git; cd dodi; pip3 install -r requirements.txt; python3 setup.py install; cd /app

# fix manager for ubuntu noble
# RUN rm /usr/lib/python3*/EXTERNALLY-MANAGED

# python depedancies
COPY requirements.txt /app
RUN pip3 install -r requirements.txt

# packages broken in noble
# RUN git clone https://github.com/arvkevi/kneed.git && cd kneed && pip3 install -e . && cd ..
# RUN git clone https://github.com/ogayot/khmer.git && cd khmer && pip3 install . && cd ..


# r dependacies
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('copynumber')"

# SPAdes
COPY SPAdes-3.15.5 /app/SPAdes-3.15.5
RUN cd /app/SPAdes-3.15.5; PREFIX=/usr/local ./spades_compile.sh; sed -i "s%#!/usr/bin/env python%#!/usr/bin/env python3%g" /usr/local/bin/spades.py

# pipelines
COPY chain_link /app/chain_link
COPY circos /app/circos
COPY cnp /app/cnp
COPY counts /app/counts
COPY RRAssembler /app/RRAssembler

# compile contig assembly
RUN cd /app/RRAssembler/; chmod +x compile_commands.sh; ./compile_commands.sh; cd /app

RUN sed -i '/disable ghostscript format types/,+6d' /etc/ImageMagick-6/policy.xml
