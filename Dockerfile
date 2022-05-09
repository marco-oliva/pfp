# Get the base Ubuntu image from Docker Hub
FROM ubuntu:22.04

# Update apps on the base image
RUN apt-get -y update && apt-get upgrade -y

# Install GCC and dependencies
RUN apt-get install -y wget bzip2 autoconf automake make cmake gcc perl zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev libssl-dev libncurses5-dev
WORKDIR /usr/src/
RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 -O htslib.tar.bz2 \
    && tar -xjvf htslib.tar.bz2 \
    && cd htslib-1.15.1 \
    && autoreconf -i \
    && ./configure \
    && make\
    && make install

# Start building
COPY . /usr/src/pfp
RUN mkdir -p /usr/src/pfp/build

# Specify the working directory
WORKDIR /usr/src/pfp/build
RUN cmake ..
RUN make -j
RUN make test
RUN make install

# Run pfp++
CMD ["pfp++ --help"]