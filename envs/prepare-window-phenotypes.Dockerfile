FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.title="mvSuSiE window phenotype preparation" \
      org.opencontainers.image.description="R environment for filtering one mvSuSiE window of gzipped phenotype BED files" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        ca-certificates \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus -1 \
    optparse \
    dplyr \
    purrr \
    readr \
    tibble

COPY scripts/prepare_trans_window.R scripts/trans_window_cli.R /opt/mvsusie/scripts/
WORKDIR /opt/mvsusie

CMD ["Rscript"]
