# Get the base Ubuntu image from Docker Hub
FROM ubuntu:22.04

# Update apps on the base image
RUN apt-get -y update && apt-get upgrade -y

# Install GCC and dependencies
RUN apt-get -y install build-essential gcc make cmake

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
CMD ["pfp++ --version"]