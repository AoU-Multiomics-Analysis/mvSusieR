FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.title="mvSuSiE trans-window fine-mapping" \
      org.opencontainers.image.description="R environment for preparing and fitting mvSuSiE in one trans window" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        ca-certificates \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        git \
        libglpk-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libssl-dev \
        libtiff5-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus -1 \
    data.table \
    optparse \
    remotes

# mvsusieR 0.3.0 requires susieR >= 0.15.54. The current susieR development
# branch is installed first because the CRAN release is older than that floor.
RUN Rscript -e 'remotes::install_github("stephenslab/susieR@master", dependencies = TRUE, upgrade = "never")' \
    && Rscript -e 'remotes::install_github("stephenslab/mvsusieR@master", dependencies = TRUE, upgrade = "never")' \
    && Rscript -e 'stopifnot(utils::packageVersion("mvsusieR") >= "0.3.0", utils::packageVersion("susieR") >= "0.15.54")'

COPY scripts/trans_window_io.R \
     scripts/trans_window_preprocess.R \
     scripts/trans_window_model.R \
     scripts/trans_window_cli.R \
     scripts/run_window_mvsusie.R \
     scripts/summarize_window.R \
     scripts/merge_window_outputs.R \
     /opt/mvsusie/scripts/

WORKDIR /opt/mvsusie

CMD ["Rscript"]
