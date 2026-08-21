FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="TransQTL BAM filtering and RNA-SeQC2" \
      org.opencontainers.image.description="samtools and RNA-SeQC2 environment for TransQTLBamFiltering" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="BSD-3-Clause"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/transqtl-bam-filtering.environment.yml /tmp/environment.yml

RUN micromamba config set channel_priority strict \
    && micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && samtools --version \
    && rnaseqc --version

WORKDIR /data

CMD ["bash"]
