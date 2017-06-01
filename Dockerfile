FROM ubuntu:14.04

RUN apt-get update
RUN apt-get upgrade -y

# install python and PIP
RUN apt-get install python2.7 python2.7-dev -y
RUN apt-get install cython -y
RUN apt-get install python-setuptools -y
RUN apt-get install python-pip -y
RUN apt-get install python-numpy -y
RUN apt-get install zlib1g-dev -y
# ideally, we don't want to use the ubuntu version
# here because this bloats the Docker image
# RUN apt-get install python-pandas -y
RUN easy_install -U distribute

# numpy and such
RUN apt-get install build-essential -y
RUN apt-get install gfortran -y
RUN apt-get install -y libatlas-base-dev
RUN apt-get install pkg-config -y
RUN apt-get install software-properties-common python-software-properties -y
RUN apt-get install cmake -y

# hap.py v0.3.7+ requires g++ 4.9
RUN add-apt-repository ppa:ubuntu-toolchain-r/test -y
RUN apt-get update
RUN apt-get install g++-4.9 -y
RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 50
RUN update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 50

# upgrade to newest versions / install
RUN pip install --upgrade cython
RUN pip install --upgrade numpy
RUN pip install --upgrade pandas
RUN pip install pybedtools
RUN pip install pysam
RUN pip install bx-python
RUN pip install scipy

# boost and cmake for hap.py
RUN apt-get install libncurses5-dev -y
RUN apt-get install git -y
RUN apt-get install samtools -y
RUN apt-get install bzip2 -y
RUN apt-get install wget -y
RUN apt-get install libbz2-dev -y


RUN echo oracle-java8-installer shared/accepted-oracle-license-v1-1 select true | debconf-set-selections && \
    add-apt-repository -y ppa:webupd8team/java && \
    apt-get update && \
    apt-get install -y oracle-java8-installer && \
    rm -rf /var/lib/apt/lists/* && \
    rm -rf /var/cache/oracle-jdk8-installer

RUN apt-get clean -y

# copy git repository into the image
RUN mkdir -p /opt/hap.py-source
COPY . /opt/hap.py-source/

# make minimal HG19 reference sequence
RUN mkdir -p /opt/hap.py-data

# download HG19 reference data
WORKDIR /opt/hap.py-data

# This downloads from UCSC
RUN /opt/hap.py-source/src/sh/make_hg19.sh

# this downloads our copy from box.com (once available)
# RUN wget https://illumina.app.box.com/s/1vwpu1w23p2qupukjjlrenq81mrtth1z/hg19.fa.bz2 -O /opt/hap.py-source/src/data/hg19.fa.bz2
# RUN bzcat /opt/hap.py-source/src/data/hg19.fa.bz2 > /opt/hap.py-data/hg19.fa
RUN samtools faidx /opt/hap.py-data/hg19.fa

# get + install ant
WORKDIR /opt
RUN wget http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.7-bin.tar.gz && \
    tar xzf apache-ant-1.9.7-bin.tar.gz && \
    rm apache-ant-1.9.7-bin.tar.gz
ENV PATH $PATH:/opt/apache-ant-1.9.7/bin

# run hap.py installer in the image
WORKDIR /opt/hap.py-source
RUN HG19=/opt/hap.py-data/hg19.fa python install.py /opt/hap.py --with-rtgtools
WORKDIR /
RUN rm -rf /opt/hap.py-source

# Optionally set up I/O folders
# RUN mkdir /data && chgrp users /data && chmod a+rwx /data
# VOLUME /data
# RUN useradd happy -m -s /bin/bash -G users
# USER happy
