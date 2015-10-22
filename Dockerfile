#############################################################
# Dockerfile based on Ubuntu Image
# supporting Decock et al. data descriptor paper
#############################################################

# Set the base image to use to Ubuntu
FROM ubuntu:14.04

# Set the file maintainer (your name - the file's author)
MAINTAINER "Maté Ongenaert" mate.ongenaert@gmail.com

RUN apt-get update && apt-get install -y build-essential \
								  zlib1g-dev \
								  zlibc \
								  openjdk-7-jre \
								  git \
								  libboost-dev \
								  autoconf \
								  libncursesw5-dev \
								  libncurses5 \
								  ncurses-dev \
								  libboost-thread-dev \
								  python3-pip \
								  samtools \
								  unzip \
									python \
									curl \
									bedtools

RUN mkdir /opt/software


#Picard, bowtie, SRAToolkit, MACS

ADD http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.4.zip /opt/software/
ADD http://downloads.sourceforge.net/project/samstat/samstat-1.5.1.tar.gz /opt/software/
ADD http://genome.sph.umich.edu/w/images/7/70/BamUtilLibStatGen.1.0.13.tgz /opt/software/
ADD https://github.com/broadinstitute/picard/releases/download/1.140/picard-tools-1.140.zip /opt/software/
ADD https://github.com/BenLangmead/bowtie2/releases/download/v2.2.6/bowtie2-2.2.6-linux-x86_64.zip /opt/software/
ADD http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.4-1/sratoolkit.2.5.4-1-ubuntu64.tar.gz /opt/software/
ADD https://pypi.python.org/packages/source/M/MACS/MACS-1.4.3.tar.gz /opt/software/


WORKDIR /opt/software

RUN unzip fastqc_v0.11.4.zip
RUN rm fastqc_v0.11.4.zip

RUN tar zxvf samstat-1.5.1.tar.gz && mv samstat-1.5.1 samstat
RUN rm samstat-1.5.1.tar.gz

RUN tar zxvf BamUtilLibStatGen.1.0.13.tgz && mv bamUtil_1.0.13 bamUtil
RUN rm BamUtilLibStatGen.1.0.13.tgz

RUN unzip bowtie2-2.2.6-linux-x86_64.zip && mv bowtie2-2.2.6 bowtie2
RUN rm bowtie2-2.2.6-linux-x86_64.zip

WORKDIR /opt/software/bowtie2/
RUN curl ftp://ftp.ccb.jhu.edu/pub/data/bowtie2_indexes/hg19.zip -o hg19.zip
RUN rm hg19.zip

WORKDIR /opt/software

RUN unzip picard-tools-1.140.zip && mv picard-tools-1.140 picard
RUN rm picard-tools-1.140.zip

RUN tar zxvf sratoolkit.2.5.4-1-ubuntu64.tar.gz && mv sratoolkit.2.5.4-1-ubuntu64 sratoolkit
RUN rm sratoolkit.2.5.4-1-ubuntu64.tar.gz

RUN tar zxvf MACS-1.4.3.tar.gz && mv MACS-1.4.3 MACS
RUN rm MACS-1.4.3.tar.gz

WORKDIR /opt/software/MACS

RUN python setup.py install --prefix /opt/software/MACS/

WORKDIR /opt/software/samstat/

RUN ./configure
RUN make

WORKDIR /opt/software/bamUtil/

RUN make all


WORKDIR /opt/software

ENV PYTHONPATH=/opt/software/MACS/lib/python2.7/site-packages/:$PYTHONPATH
ENV PATH=/opt/software/bamUtil/bamUtil/bin/:/opt/software/samstat/src/:/opt/software/FastQC:/opt/software/bowtie2:/opt/software/sratoolkit/bin/:/opt/software/picard/:/opt/software/MACS/bin/:$PATH

ENV LANG en_US.UTF-8
# Default action
CMD ["/bin/bash"]
