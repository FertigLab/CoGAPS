# --platform=linux/amd64 to avoid 'no match for platform in the manifest' on M1
FROM rocker/tidyverse:4

COPY . /cogaps
WORKDIR /cogaps

RUN sudo apt-get update -y && \
    apt-get upgrade -y && \
    apt-get install libhdf5-dev build-essential patch cmake -y

RUN Rscript -e 'devtools::install_deps()'

#https://github.com/r-lib/devtools/issues/2395
RUN Rscript -e 'devtools::install()'