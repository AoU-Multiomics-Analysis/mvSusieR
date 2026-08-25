FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="mvSuSiE trans-window fine-mapping" \
      org.opencontainers.image.description="R environment for preparing and fitting mvSuSiE in one trans window" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/trans-window-mvsusie.environment.yml /tmp/environment.yml

# Activate the base environment for image-build commands and task runtimes.
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN micromamba config set channel_priority strict \
    && micromamba install --yes --name base --override-channels --strict-channel-priority --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && Rscript -e 'stopifnot(packageVersion("mvsusieR") >= package_version("0.3.0"), packageVersion("susieR") >= package_version("0.15.0"), requireNamespace("mashr", quietly = TRUE), is.function(mvsusieR::mvsusie_plot))'

COPY --chown=$MAMBA_USER:$MAMBA_USER scripts/trans_window_io.R \
     scripts/trans_window_logging.R \
     scripts/trans_window_preprocess.R \
     scripts/trans_window_model.R \
     scripts/trans_window_prior.R \
     scripts/trans_window_cli.R \
     scripts/prepare_window.R \
     scripts/fit_window.R \
     scripts/run_window_mvsusie.R \
     scripts/summarize_window.R \
     scripts/merge_window_outputs.R \
     scripts/plot_window_mvsusie.R \
     /opt/mvsusie/scripts/

WORKDIR /opt/mvsusie

CMD ["Rscript"]
