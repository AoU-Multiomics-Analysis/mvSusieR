FROM ubuntu:24.04

LABEL org.opencontainers.image.title="mvSuSiE window genotype preparation" \
      org.opencontainers.image.description="Tabix environment for extracting one mvSuSiE window of dosage data" \
      org.opencontainers.image.source="https://github.com/AoU-Multiomics-Analysis/mvSusieR" \
      org.opencontainers.image.licenses="GPL-3.0-or-later"

RUN apt-get update \
    && DEBIAN_FRONTEND=noninteractive apt-get install --yes --no-install-recommends \
        ca-certificates \
        tabix \
    && rm -rf /var/lib/apt/lists/*

CMD ["tabix"]
