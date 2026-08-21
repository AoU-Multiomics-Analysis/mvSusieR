FROM ubuntu:24.04

ARG RNASEQC_VERSION=2.4.2

LABEL org.opencontainers.image.title="TransQTL BAM filtering and RNA-SeQC2" \
      org.opencontainers.image.description="samtools and RNA-SeQC2 environment for TransQTLBamFiltering" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="BSD-3-Clause"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        ca-certificates \
        curl \
        gzip \
        samtools \
    && curl --fail --location --retry 3 \
        "https://github.com/getzlab/rnaseqc/releases/download/v${RNASEQC_VERSION}/rnaseqc.v${RNASEQC_VERSION}.linux.gz" \
        --output "/tmp/rnaseqc.v${RNASEQC_VERSION}.linux.gz" \
    && gzip --decompress "/tmp/rnaseqc.v${RNASEQC_VERSION}.linux.gz" \
    && install -m 0755 "/tmp/rnaseqc.v${RNASEQC_VERSION}.linux" /usr/local/bin/rnaseqc \
    && rm -rf /var/lib/apt/lists/* "/tmp/rnaseqc.v${RNASEQC_VERSION}.linux" \
    && samtools --version \
    && rnaseqc --version

WORKDIR /data

CMD ["bash"]
