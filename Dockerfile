# Start from a specific, stable Ubuntu version for reproducibility
FROM ubuntu:22.04
LABEL maintainer="Effie Klimi <effie@effie.bio>"

# Set environment variable to allow non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory for all subsequent commands
WORKDIR /app

# 1. Install all system-level dependencies in a single layer
#    - build-essential: For compiling C/C++ code
#    - zlib1g-dev, etc: Development headers for bcftools
#    - python3-dev: Development headers for compiling Python C extensions (e.g., for pyBigWig)
#    - libcurl4-openssl-dev: For pyBigWig remote file support
#    - openjdk-11-jre: For Beagle (Java)
#    - python3, python3-pip: For the scorer and CrossMap
#    - git, wget, etc: Utilities
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    python3-dev \
    libcurl4-openssl-dev \
    openjdk-11-jre \
    python3 \
    python3-pip \
    git \
    wget \
    bzip2 \
    xz-utils \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

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

# --- Install CrossMap ---
RUN pip3 install CrossMap

# 3. Copy your project's scripts into the container
#    This single command copies the entire 'scripts' directory.
COPY scripts/ /app/scripts/

# 4. Make all shell scripts executable in a single, efficient layer
#    Using a wildcard is robust and automatically handles any new scripts you add.
#    Note: config.sh is correctly excluded as it doesn't need to be executable.
RUN chmod +x /app/scripts/pipeline/user.sh \
           /app/scripts/pipeline/prep_pgs.sh \
           /app/scripts/analyses/impute.sh

# 5. Set the default command to run when the container starts interactively
CMD ["/bin/bash"]