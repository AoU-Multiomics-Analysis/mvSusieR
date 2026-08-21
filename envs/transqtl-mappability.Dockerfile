FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="TransQTL mappability generation" \
      org.opencontainers.image.description="GenMap and BED utilities for read-length-specific TransQTL mappability masks" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="BSD-3-Clause"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/transqtl-mappability.environment.yml /tmp/environment.yml

# Set up the base environment for both image build commands and task runtimes.
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN micromamba config set channel_priority strict \
    && micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && genmap --version \
    && bedtools --version

WORKDIR /data

CMD ["bash"]
