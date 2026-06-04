FROM rocker/tidyverse:4

COPY . /cogaps
WORKDIR /cogaps

RUN sudo apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install -y libhdf5-dev build-essential patch libuv1 libuv1-dev

RUN Rscript -e 'install.packages(c("remotes","BiocManager","devtools"), repos="https://cloud.r-project.org")'

RUN Rscript -e 'BiocManager::install(c("S4Vectors", "SingleCellExperiment", "SummarizedExperiment", "rhdf5", "fgsea", "sparseMatrixStats"), ask=FALSE)'

RUN Rscript -e 'remotes::install_deps(".", dependencies=TRUE)'

RUN Rscript -e 'devtools::install(".", dependencies=TRUE)'