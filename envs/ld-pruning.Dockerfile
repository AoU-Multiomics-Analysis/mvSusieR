FROM ubuntu:24.04

ARG PLINK2_DATE=20260818

LABEL org.opencontainers.image.title="mvSusieR LD pruning" \
      org.opencontainers.image.description="PLINK 2 container for the mvSusieR LD-pruning workflow" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        ca-certificates \
        curl \
        unzip \
    && curl --fail --location --retry 3 \
        "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_${PLINK2_DATE}.zip" \
        --output /tmp/plink2.zip \
    && unzip -q /tmp/plink2.zip -d /tmp/plink2 \
    && install -m 0755 /tmp/plink2/plink2 /usr/local/bin/plink2 \
    && rm -rf /var/lib/apt/lists/* /tmp/plink2 /tmp/plink2.zip

CMD ["plink2"]
