#Download base image ubuntu
FROM combinelab/salmon

# Install needed libraries
ENV PACKAGES git gcc make g++ libboost-all-dev liblzma-dev libbz2-dev \
    ca-certificates zlib1g-dev libcurl4-openssl-dev curl unzip autoconf apt-transport-https ca-certificates gnupg software-properties-common wget default-jdk python3-pip cmake vim
ENV SALMON_VERSION 1.10.0
ENV DST=/tmp
ENV URL=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
ENV ZIP=fastqc_v0.11.5.zip


# Update Ubuntu Software repository
RUN apt-get update && \
    apt-get install -y ${PACKAGES} && \
    apt-get clean

# FASTQC
RUN wget $URL/$ZIP -O $DST/$ZIP && \
  unzip - $DST/$ZIP -d $DST && \
  rm $DST/$ZIP && \
  cd $DST/FastQC && \
  chmod 755 fastqc && \
  ln -s $DST/FastQC/fastqc /usr/local/bin/fastqc
  
# MultiQC
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install --upgrade Pillow
RUN pip3 install multiqc

ENV PATH /usr/local/bin:$PATH
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8

WORKDIR /data
