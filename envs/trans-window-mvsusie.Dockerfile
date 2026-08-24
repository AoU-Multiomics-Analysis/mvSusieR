FROM rocker/r-ver:4.4.1

LABEL org.opencontainers.image.title="mvSuSiE trans-window fine-mapping" \
      org.opencontainers.image.description="R environment for preparing and fitting mvSuSiE in one trans window" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        ca-certificates \
        build-essential \
        gfortran \
        libcurl4-openssl-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        git \
        libglpk-dev \
        libgsl-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        ripgrep \
        libssl-dev \
        libtiff5-dev \
        libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN install2.r --error --skipinstalled --ncpus -1 \
    data.table \
    dplyr \
    ggplot2 \
    optparse \
    purrr \
    readr \
    remotes \
    R.utils \
    stringr \
    tibble

# Install exact source revisions. Do not install suggested documentation and test packages.
RUN Rscript -e 'remotes::install_github("stephenslab/susieR@65f3586a865fb6748cb4f9df50510ac577706348", dependencies = c("Depends", "Imports", "LinkingTo"), upgrade = "never")' \
    && Rscript -e 'install.packages("mashr", repos = "https://cloud.r-project.org")' \
    && Rscript -e 'remotes::install_github("stephenslab/mvsusieR@ebd1133953005fa70c6b338727b5fe9222e2a1c2", dependencies = c("Depends", "Imports", "LinkingTo"), upgrade = "never")' \
    && Rscript -e 'stopifnot(requireNamespace("ggplot2", quietly = TRUE), requireNamespace("mashr", quietly = TRUE), packageDescription("susieR")$RemoteSha == "65f3586a865fb6748cb4f9df50510ac577706348", packageDescription("mvsusieR")$RemoteSha == "ebd1133953005fa70c6b338727b5fe9222e2a1c2")'

COPY scripts/trans_window_io.R \
     scripts/trans_window_logging.R \
     scripts/trans_window_preprocess.R \
     scripts/trans_window_model.R \
     scripts/trans_window_prior.R \
     scripts/trans_window_cli.R \
     scripts/fit_window.R \
     scripts/run_window_mvsusie.R \
     scripts/summarize_window.R \
     scripts/merge_window_outputs.R \
     scripts/plot_window_mvsusie.R \
     /opt/mvsusie/scripts/

WORKDIR /opt/mvsusie

CMD ["Rscript"]
