FROM ubuntu:18.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y \
        autoconf \
        build-essential \
        bzip2 \
        cmake \
        cython \
        git \
        libbz2-dev \
        libncurses5-dev \
        openjdk-8-jdk \
        pkg-config \                    
        python \
        python2.7 \
        python2.7-dev \                    
        python-setuptools \
        python-pip \
        python-psutil \                    
        python-numpy \
        python-pandas \
        python-distribute \
        python-pysam \
        python-scipy \                    
        software-properties-common \
        wget \
        zlib1g-dev && \
    apt-get clean -y

RUN pip install bx-python

# copy git repository into the image
RUN mkdir -p /opt/hap.py-source
COPY . /opt/hap.py-source/

# make minimal HG19 reference sequence
RUN mkdir -p /opt/hap.py-data

# download HG19 reference data
WORKDIR /opt/hap.py-data

# get + install ant
WORKDIR /opt
RUN wget http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.7-bin.tar.gz && \
    tar xzf apache-ant-1.9.7-bin.tar.gz && \
    rm apache-ant-1.9.7-bin.tar.gz
ENV PATH $PATH:/opt/apache-ant-1.9.7/bin

# run hap.py installer in the image, don't run tests since we don't have a reference file
WORKDIR /opt/hap.py-source
RUN python install.py /opt/hap.py --with-rtgtools --no-tests
WORKDIR /opt/hap.py

# run basic tests
RUN bin/test_haplotypes

# remove source folder
WORKDIR /
RUN rm -rf /opt/hap.py-source

