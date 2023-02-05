# Get the base Ubuntu image from Docker Hub
FROM ubuntu:20.04

# Install GCC and dependencies
RUN apt-get -y update \
    && DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata \
    && apt-get install -y git wget bzip2 autoconf automake make cmake gcc g++ perl zlib1g-dev libbz2-dev liblzma-dev \
    libcurl4-gnutls-dev libssl-dev libncurses5-dev

WORKDIR /usr/src/
RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 -O htslib.tar.bz2 \
    && tar -xjvf htslib.tar.bz2 \
    && cd htslib-1.15.1 \
    && autoreconf -i \
    && ./configure \
    && make \
    && make install

# Start building
COPY . /usr/src/pfp
WORKDIR /usr/src/pfp/build
RUN cmake -DENABLE_MIMALLOC=ON .. \
    && make -j

# Get binaries
WORKDIR /pfp/bin
RUN cp \
    /usr/src/pfp/build/pfp++ \
    /usr/src/pfp/build/pfp++64 \
    /usr/src/pfp/build/vcf_to_fa \
    /usr/src/pfp/build/mpfp++ \
    /usr/src/pfp/build/mpfp++64 \
    /usr/src/pfp/build/check \
    /usr/src/pfp/build/check64 \
    /usr/src/pfp/build/exprop \
    /usr/src/pfp/build/exprop64  \
    .
ENV PATH /pfp/bin:$PATH