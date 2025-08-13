FROM --platform=linux/amd64 ubuntu:22.04
LABEL maintainer="Effie Klimi <effie@effie.bio>"

# Set environment variable to allow non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory for all subsequent commands
WORKDIR /app

# 1. Install all system-level dependencies in a single layer
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python3-dev \
    libcurl4-openssl-dev \
    openjdk-17-jre \
    python3 \
    python3-pip \
    git \
    wget \
    bzip2 \
    unzip \
    xz-utils \
    ca-certificates \
    samtools \
    tabix \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

RUN ln -s /usr/bin/python3 /usr/bin/python

# 2. Install key bioinformatics tools by downloading binaries/source

# --- Install bcftools & tabix by compiling from source ---
RUN wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2 \
    && tar -xjvf bcftools-1.18.tar.bz2 \
    && cd bcftools-1.18 \
    && ./configure && make && make install \
    && cd .. \
    && rm -r bcftools-1.18 bcftools-1.18.tar.bz2

# --- Install Eagle (with proper permissions and verification) ---
RUN wget https://data.broadinstitute.org/alkesgroup/Eagle/downloads/Eagle_v2.4.1.tar.gz \
    && tar -zxvf Eagle_v2.4.1.tar.gz \
    && ls -la Eagle_v2.4.1/ \
    && chmod +x Eagle_v2.4.1/eagle \
    && mv Eagle_v2.4.1/eagle /usr/local/bin/ \
    && rm -r Eagle_v2.4.1.tar.gz Eagle_v2.4.1 \
    && which eagle \
    && eagle --help || echo "Eagle installed but help command failed (this might be normal)"

# --- Install CrossMap (improved version with proper PATH handling) ---
RUN pip3 install --no-cache-dir CrossMap \
    && echo "=== Debugging CrossMap installation ===" \
    && find /usr -name "*CrossMap*" -type f 2>/dev/null || echo "No CrossMap files found in /usr" \
    && find /home -name "*CrossMap*" -type f 2>/dev/null || echo "No CrossMap files found in /home" \
    && find /root -name "*CrossMap*" -type f 2>/dev/null || echo "No CrossMap files found in /root" \
    && ls -la /usr/local/bin/ | grep -i cross || echo "No CrossMap in /usr/local/bin" \
    && ls -la ~/.local/bin/ 2>/dev/null | grep -i cross || echo "No ~/.local/bin or no CrossMap there" \
    && python3 -c "import CrossMap; print('CrossMap Python module installed successfully')" \
    && python3 -c "import CrossMap; print(CrossMap.__file__)" \
    && echo "=== Creating CrossMap wrapper ===" \
    && echo '#!/bin/bash\npython3 -m CrossMap "$@"' > /usr/local/bin/CrossMap.py \
    && chmod +x /usr/local/bin/CrossMap.py \
    && echo "=== Testing CrossMap ===" \
    && CrossMap.py --help || python3 -m CrossMap --help || echo "CrossMap help failed but may still work"

# --- Install PLINK 2 (Linux x86_64 static binary) ---
# Pin a known-good build; update URL if you want a newer build later.
ENV PLINK2_ZIP=plink2_linux_x86_64_20250806.zip
RUN wget -q https://s3.amazonaws.com/plink2-assets/alpha6/${PLINK2_ZIP} \
    && unzip -q ${PLINK2_ZIP} \
    && mv plink2 /usr/local/bin/ \
    && chmod +x /usr/local/bin/plink2 \
    && rm -f ${PLINK2_ZIP} plink2.license plink2.md \
    && /usr/local/bin/plink2 --version || /usr/local/bin/plink2 -v || echo "plink2 installed (version check non-fatal)"

# --- Python packages (single layer, pinned where sensible) ---
RUN pip3 install --no-cache-dir \
    numpy==1.26.4 \
    pandas==2.2.2 \
    requests==2.32.3 \
    pyfaidx==0.8.1.1 \
    pyliftover==0.4 \
    scipy \
    scikit-learn \
    cyvcf2 \
    joblib \
 && python3 - <<'PY'
import numpy, pandas, requests
from pyfaidx import Fasta
try:
    from pyliftover import LiftOver
except Exception as e:
    raise SystemExit("pyliftover import failed: %s" % e)
print("Imports OK")
PY


# 3. Copy your project's scripts into the container
COPY scripts/ /app/scripts/

# 4. Make all shell scripts executable in a single, efficient layer
RUN chmod +x /app/scripts/pipeline/user.sh \
           /app/scripts/pipeline/add_pgs.sh \
           /app/scripts/pipeline/setup.sh \
           /app/scripts/helpers/pruned_panel.sh \
           /app/scripts/pipeline/qc_genome.sh \
           /app/scripts/helpers/download_data.sh \
           /app/scripts/helpers/build_refs.sh \
           /app/scripts/helpers/vcf_to_pfile.sh \
           /app/scripts/helpers/per_ancestry_maf.sh \
           /app/scripts/analyses/impute.sh

# 5. Set the default command to run when the container starts interactively
CMD ["/bin/bash"]
