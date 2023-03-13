FROM gitlab-registry.internal.sanger.ac.uk/sanger-pathogens/docker-images/pathogens-base:0.2
WORKDIR /opt

ARG DEBIAN_FRONTEND=noninteractive
ARG MULTIMAPPER_BUILD_DIR=multimapper

RUN apt update -qq -y && \
    apt install -y \
    python3 \
    python3-setuptools \
    python3-pip \
    wget \
    make \
    git \
    g++ \
    meson \
    samtools \
    bcftools \
    freebayes \
    flex \
    bison \
    libgmp3-dev \
    cmake \
    libgtest-dev \
    && apt clean -y

#Install smalt
RUN wget http://sourceforge.net/projects/smalt/files/smalt-0.7.6-static.tar.gz/download \
        && tar -zxvf download \
        && cd smalt-0.7.6 \
        && ./configure \
        && make \
        && make install \
        && cd ~
ENV PATH="/opt/smalt-0.7.6:$PATH"

#Install BWA-mem2
RUN git clone https://github.com/bwa-mem2/bwa-mem2 \
        && cd bwa-mem2 \
        && git submodule init \
        && git submodule update \
        && make \
        && cd ~
ENV PATH="$PATH:/opt/bwa-mem2/"

#Install RAxML-NG
RUN wget https://github.com/amkozlov/raxml-ng/archive/refs/tags/1.1.0.tar.gz \
        && tar -xvf 1.1.0.tar.gz \
        && cd raxml-ng-1.1.0 \
        && mkdir build \
        && cd build \
        && cmake .. \
        && make \
        && cd ~
# Put binary on PATH
ENV PATH=/opt/raxml-ng:${PATH}

#Install java-19
RUN wget https://download.oracle.com/java/19/latest/jdk-19_macos-x64_bin.tar.gz \
        && tar -xvf ~/Downloads/openjdk-16.0.1_linux-x64_bin.tar.gz -C /opt \
        && update-alternatives --install /usr/bin/java java /opt/jdk-16.0.1/bin/java 1000 \
        && update-alternatives --install /usr/bin/javac javac /opt/jdk-16.0.1/bin/javac 1000 \
        && update-alternatives --config java \
        && update-alternatives --config javac

#Install Sambamba
#ARG LDC=1.30.0
#RUN wget https://github.com/ldc-developers/ldc/releases/download/v$LDC/ldc2-1.7.0-linux-x86_64.tar.xz \
#        && tar xvJf ldc2-1.7.0-linux-x86_64.tar.xz \
#        && export PATH=$HOME/ldc2-1.7.0-linux-x86_64/bin:$PATH \
#        && export LIBRARY_PATH=$HOME/ldc2-1.7.0-linux-x86_64/lib \
#        && git clone --recursive https://github.com/biod/sambamba.git \
#        && cd sambamba \
#        && make

#Install Bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.5.1.tar.gz \
    && tar -xvf v2.5.1.tar.gz \
    && cd bowtie2-v2.5.1 \
    && make \
    && make static-libs \
    && make STATIC_BUILD=1
ENV PATH=/opt/bowtie2-v2.5.1

# Install multimapper
RUN mkdir -p $MULTIMAPPER_BUILD_DIR
COPY . $MULTIMAPPER_BUILD_DIR
RUN cd $MULTIMAPPER_BUILD_DIR \
    && python3 setup.py test \
    && python3 setup.py install