FROM rocker/tidyverse:4

COPY . /cogaps
WORKDIR /cogaps

RUN sudo apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y \
      libhdf5-dev build-essential patch \
      libuv1 libuv1-dev \
      autoconf automake libtool autoconf-archive

RUN autoreconf -fi

RUN Rscript -e 'install.packages("pak")' && \
    Rscript -e 'pak::local_install_deps(".")' && \
    Rscript -e 'pak::pkg_install("sparseMatrixStats")' # sparseMatrixStats is for the nextflow module

RUN R CMD INSTALL .
