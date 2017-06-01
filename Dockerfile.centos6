FROM centos:6.7

ENV HOME /root
ENV LANG en_US.UTF-8
ENV LC_ALL en_US.UTF-8

RUN  yum install -y wget curl tar gzip bzip2 which perl

RUN  mkdir -p /opt/illumina/haplocompare && \
     cd /opt/illumina/haplocompare && \
     wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
     bash Miniconda2-latest-Linux-x86_64.sh -b -p /opt/miniconda2

ENV PATH /opt/miniconda2/bin:${PATH}

RUN  yum install -y centos-release-scl && \
     yum install -y devtoolset-3-gcc devtoolset-3-binutils && \
     yum install -y devtoolset-3-gcc-c++ devtoolset-3-gcc-gfortran && \
     yum install -y cmake bzip2-devel ncurses-devel zlib-devel git

RUN  conda config --add channels bioconda && \
     conda install -y cython numpy scipy biopython matplotlib pandas pysam bx-python pyvcf cyvcf2 nose

RUN yum upgrade -y && yum update -y && yum clean all

COPY . /opt/illumina/haplocompare/hap.py-source
WORKDIR /opt/illumina/haplocompare/hap.py-source

# patch samtools for centos 6
RUN  cd external && \
     tar xvzf samtools.tar.gz && \
     cd samtools && \
     cat Makefile | sed 's/-ldl/-ldl -ltinfo/' > Makefile.bak && \
     mv -f Makefile.bak Makefile && \
     cd .. && \
     rm -f samtools.tar.gz && \
     tar czvf samtools.tar.gz samtools


ENV JDK_VERSION 8u11
ENV JDK_BUILD_VERSION b12
RUN curl -LO "http://download.oracle.com/otn-pub/java/jdk/$JDK_VERSION-$JDK_BUILD_VERSION/jdk-$JDK_VERSION-linux-x64.rpm" -H 'Cookie: oraclelicense=accept-securebackup-cookie' && rpm -i jdk-$JDK_VERSION-linux-x64.rpm && rm -f jdk-$JDK_VERSION-linux-x64.rpm && yum clean all
ENV JAVA_HOME /usr/java/default

# get + install ant
WORKDIR /opt
RUN wget http://archive.apache.org/dist/ant/binaries/apache-ant-1.9.7-bin.tar.gz && \
    tar xzf apache-ant-1.9.7-bin.tar.gz && \
    rm apache-ant-1.9.7-bin.tar.gz
ENV PATH $PATH:/opt/apache-ant-1.9.7/bin

RUN  mkdir -p /opt/illumina/haplocompare/data && \
     cd /opt/illumina/haplocompare/data && \
     bash /opt/illumina/haplocompare/hap.py-source/src/sh/make_hg19.sh

ENV HG19 /opt/illumina/haplocompare/data/hg19.fa

WORKDIR /opt/illumina/haplocompare/hap.py-source
RUN scl enable devtoolset-3 'python install.py /opt/illumina/haplocompare/hap.py --with-rtgtools'
