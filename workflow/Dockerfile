# Dockerfile for UK-Essen COVID Sequencing pipeline

FROM ubuntu:20.10
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    git \
    make \
    libgsl-dev\
    build-essential\
    libz-dev \
    bwa \
    samtools \
    spades \
    autotools-dev \
    automake \
    picard-tools \
    unzip \
    default-jre \
    wget \
    cmake \
    curl \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev\
    libssl-dev


LABEL maintainer=folker.meyer@uk-essen.de

#### install BWE-MEM
RUN cd /root  \
    && git clone --recursive https://github.com/bwa-mem2/bwa-mem2 \
    && cd bwa-mem2 && make \
    && install -m755 -s ./bwa-mem2 /usr/local/bininstall \
    && cd /root ; rm -rf bwa-mem2


#### install fastp
RUN cd /root \
    && git clone https://github.com/OpenGene/fastp.git \
    && cd fastp \
    && make \
    && make install \
    && cd /root \
    && rm -rf fastp

#### install minimap2
RUN cd /root \
    && git clone https://github.com/lh3/minimap2 \
    && cd minimap2 \
    && make \
    #&& install -m755 -s ./minimapbwa-mem2 /usr/local/bin \
    && cd /root \
    && rm -rf minimap2

#### install HTSLIB required by lofreq
RUN cd /root \
    && wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2  \
    && tar -vxjf htslib-1.9.tar.bz2 \
    && rm -f htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && autoreconf \
    && make  \
    && make install

#### install lofreq
RUN cd /root \
    && git clone https://github.com/CSB5/lofreq \
    && cd lofreq \
    && ./bootstrap \
    && ./configure --with-htslib=/usr/local/lib/ \
    && make install \
    && cd /root \
    && rm -rf htslib-1.9 lofreq

#### install SnpEff
RUN cd /root \
    && wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip \
    && unzip snpEff_latest_core.zip \
    && cd snpEff \
    && echo "snpEff is install in /root/snpEff run via java -jar snpEff.jar"

#### install fastv
RUN cd /root \
    && git clone https://github.com/OpenGene/fastv \
    && cd fastv \
    && make \
    && make install \
    && cd /root \
    && rm -rf fastv

#### install megahit
RUN cd /root \
    &&  git clone https://github.com/voutcn/megahit.git \
    && cd megahit \
    && git submodule update --init \
    && mkdir build && cd build \
    && cmake .. -DCMAKE_BUILD_TYPE=Release  -DCMAKE_INSTALL_PREFIX=/usr/local \
    && make -j4
      
