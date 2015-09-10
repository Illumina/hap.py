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

# upgrade to newest versions / install
RUN pip install --upgrade cython
RUN pip install --upgrade numpy
RUN pip install --upgrade pandas
RUN pip install pybedtools
RUN pip install pysam
RUN pip install bx-python

# boost and cmake for hap.py
RUN apt-get install libncurses5-dev -y
RUN apt-get install git -y
RUN apt-get install samtools -y
RUN apt-get install bzip2 -y
RUN apt-get install wget -y
RUN apt-get install libbz2-dev -y
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

# run hap.py installer in the image
WORKDIR /opt/hap.py-source
RUN HG19=/opt/hap.py-data/hg19.fa python install.py /opt/hap.py
WORKDIR /
RUN rm -rf /opt/hap.py-source

# Optionally set up I/O folders
# RUN mkdir /data && chgrp users /data && chmod a+rwx /data
# VOLUME /data
# RUN useradd happy -m -s /bin/bash -G users
# USER happy
