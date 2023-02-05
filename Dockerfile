# Get the base Ubuntu image from Docker Hub
FROM ubuntu:20.04

# Update apps on the base image
RUN apt-get -y update

# Install GCC and dependencies
RUN DEBIAN_FRONTEND=noninteractive TZ=Etc/UTC apt-get -y install tzdata
RUN apt-get update \
    && apt-get install -y wget bzip2 autoconf automake make cmake gcc g++ perl zlib1g-dev libbz2-dev liblzma-dev \
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
RUN cmake .. \
    && make -j

# Get binaries
WORKDIR /pfp/bin
RUN cp /usr/src/pfp/build/pfp++ /usr/src/pfp/build/pfp++64 .
ENV PATH /pfp/bin:$PATH