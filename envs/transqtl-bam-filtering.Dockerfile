FROM mambaorg/micromamba:2.3.3

LABEL org.opencontainers.image.title="TransQTL BAM filtering and RNA-SeQC2" \
      org.opencontainers.image.description="samtools and RNA-SeQC2 environment for TransQTLBamFiltering" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="BSD-3-Clause"

COPY --chown=$MAMBA_USER:$MAMBA_USER envs/transqtl-bam-filtering.environment.yml /tmp/environment.yml

# Set up the base environment for both image build commands and task runtimes.
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin

RUN micromamba config set channel_priority strict \
    && micromamba install --yes --name base --file /tmp/environment.yml \
    && micromamba clean --all --yes \
    && samtools --version \
    && rnaseqc --version

WORKDIR /data

CMD ["bash"]
