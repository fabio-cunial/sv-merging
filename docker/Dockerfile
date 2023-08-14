FROM continuumio/miniconda3
MAINTAINER Fabio Cunial
ARG work_dir=/sv-merging
WORKDIR ${work_dir}


# --------------------------------- Versions -----------------------------------
ARG gcloud_version=405.0.1
ARG htslib_version=1.18
ARG samtools_version=1.18
ARG bcftools_version=1.18
ARG truvari_version=3.5.0
ARG survivor_version=1.0.7
ARG jasmine_version=1.1.5
ARG svtools_version=0.5.1
ARG svimmer_version=0.1
# ------------------------------------------------------------------------------


# OS
RUN apt-get -qqy update --fix-missing \
    && apt-get -qqy dist-upgrade \
    && apt-get install -y --no-install-recommends \
        zlib1g-dev \
        liblzma-dev \
        libbz2-dev \
        libdeflate-dev \
        libreadline-dev \
        libsqlite3-dev \
        libssl-dev \
        libcurl4-openssl-dev \
        libncurses5-dev \
        libncursesw5-dev \
        libffi-dev \
        liblzma-dev \
        libopenblas-dev \
        apt-transport-https \
        ca-certificates \
        tree \
        gnupg \
        ssh \
        time \
        curl \
        wget \
        autotools-dev \
        autoconf \
        automake \
        make \
        cmake \
        gcc \
        g++ \
        gfortran \
        build-essential \
        default-jdk \
        default-jre \
        git \
        bc \
        python3-pip \
        xz-utils \
        tk-dev
        

# GSUTIL
RUN pip3 uninstall -y crcmod && pip3 install --no-cache-dir -U crcmod
RUN wget https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && tar -xf google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && rm -f google-cloud-cli-${gcloud_version}-linux-x86_64.tar.gz \
    && yes | ./google-cloud-sdk/install.sh
ENV PATH=${work_dir}/google-cloud-sdk/bin:${PATH}

# HTSLIB
RUN wget https://github.com/samtools/htslib/releases/download/${htslib_version}/htslib-${htslib_version}.tar.bz2 \
    && tar xjf htslib-${htslib_version}.tar.bz2 \
    && rm htslib-${htslib_version}.tar.bz2 \
    && cd htslib-${htslib_version} \
    && ./configure \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf htslib-${htslib_version} \
    && bgzip --help

# SAMTOOLS
RUN wget https://github.com/samtools/samtools/releases/download/${samtools_version}/samtools-${samtools_version}.tar.bz2 \
    && tar xjf samtools-${samtools_version}.tar.bz2 \
    && rm samtools-${samtools_version}.tar.bz2 \
    && cd samtools-${samtools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf samtools-${samtools_version} \
    && samtools --help

# BCFTOOLS
RUN wget https://github.com/samtools/bcftools/releases/download/${bcftools_version}/bcftools-${bcftools_version}.tar.bz2 \
    && tar xjf bcftools-${bcftools_version}.tar.bz2 \
    && rm bcftools-${bcftools_version}.tar.bz2 \
    && cd bcftools-${bcftools_version} \
    && ./configure --without-curses \
    && make \
    && make install \
    && cd ${work_dir} \
    && rm -rf bcftools-${bcftools_version} \
    && bcftools --help

# TRUVARI
RUN pip install truvari==${truvari_version} \
    && truvari --help

# SURVIVOR
RUN wget https://github.com/fritzsedlazeck/SURVIVOR/archive/refs/tags/${survivor_version}.tar.gz -O SURVIVOR.tar.gz \
    && tar xzvf SURVIVOR.tar.gz \
    && rm SURVIVOR.tar.gz \
    && cd SURVIVOR-${survivor_version}/Debug/ \
    && make
ENV PATH=${work_dir}/SURVIVOR-${survivor_version}/Debug:${PATH}
RUN SURVIVOR

# SVPOP
COPY ./svpop_environment.yml .
RUN conda env create -f ./svpop_environment.yml \
    && conda clean --packages --tempfiles --logfiles
SHELL ["conda", "run", "-n", "svpop", "/bin/bash", "-c"]
RUN echo "conda init bash" > ~/.bashrc
RUN git clone --recursive https://github.com/EichlerLab/svpop.git

# JASMINE
COPY ./jasmine_environment.yml .
RUN conda env create -f ./jasmine_environment.yml \
    && conda clean --packages --tempfiles --logfiles
SHELL ["conda", "run", "-n", "jasmine", "/bin/bash", "-c"]
# Adding #!/bin/bash to the original shell script
COPY ./jasmine /opt/conda/envs/jasmine/bin/jasmine

# # PYENV -- SVTOOLS needed 2.7.13
# RUN git clone https://github.com/pyenv/pyenv.git ~/.pyenv \
#     && echo 'export PYENV_ROOT="$HOME/.pyenv"' >> ~/.bashrc \
#     && echo 'export PATH="$PYENV_ROOT/bin:$PATH"' >> ~/.bashrc \
#     && echo -e 'if command -v pyenv 1>/dev/null 2>&1; then\n eval "$(pyenv init -)"\nfi' >> ~/.bashrc
# RUN $HOME/.pyenv/bin/pyenv install 3.6.0 \
#     && $HOME/.pyenv/bin/pyenv versions

# SVIMMER
RUN pip install pysam \
    && wget https://github.com/DecodeGenetics/svimmer/archive/refs/tags/v${svimmer_version}.tar.gz -O svimmer.tar.gz \
    && tar xzvf svimmer.tar.gz \
    && rm svimmer.tar.gz \
    && chmod +x ./svimmer-${svimmer_version}/svimmer
# Fixing two bugs in SVIMMER
COPY ./svimmer ${work_dir}/svimmer-${svimmer_version}/svimmer
ENV PATH=${work_dir}/svimmer-${svimmer_version}:${PATH}
RUN svimmer -h

# # SVTOOLS -- I'm unable to install it.
# ARG svtools_version=0.5.1
# RUN $HOME/.pyenv/bin/pyenv global 2.7.13 \
#     && wget https://github.com/hall-lab/svtools/archive/refs/tags/v${svtools_version}.tar.gz -O svtools.tar.gz \
#     && tar xzvf svtools.tar.gz \
#     && rm svtools.tar.gz
# COPY setup.py ${work_dir}/svtools-${svtools_version}
# RUN cd svtools-${svtools_version} \
#     && pip install . \
#     && svtools --version \
#     && $HOME/.pyenv/bin/pyenv global system

# ADDITIONAL FILES
COPY ./preprocess_vcf.py ${work_dir}/
COPY *.java ${work_dir}/
RUN javac ${work_dir}/*.java

# BEDTOOLS
ARG bedtools_version=2.31.0
RUN wget https://github.com/arq5x/bedtools2/releases/download/v${bedtools_version}/bedtools.static \
    && mv ./bedtools.static ./bedtools \
    && chmod +x ./bedtools
