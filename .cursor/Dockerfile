# Cursor Cloud base image: micromamba + phc-meep (Meep, numpy, scipy, ...)
# Do not COPY the repo — Cursor checks out the workspace at runtime.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    curl \
    git \
    sudo \
    && rm -rf /var/lib/apt/lists/*

RUN id -u ubuntu &>/dev/null || useradd -m -s /bin/bash ubuntu
RUN echo "ubuntu ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers.d/ubuntu

USER ubuntu
WORKDIR /home/ubuntu

ENV MAMBA_ROOT_PREFIX=/home/ubuntu/.micromamba
ENV PATH="/home/ubuntu/.local/bin:${PATH}"

RUN mkdir -p /home/ubuntu/.local/bin \
    && curl -Ls "https://micro.mamba.pm/api/micromamba/linux-64/latest" \
      | tar -xj -C /tmp \
    && mv /tmp/bin/micromamba /home/ubuntu/.local/bin/micromamba \
    && chmod +x /home/ubuntu/.local/bin/micromamba

RUN eval "$(/home/ubuntu/.local/bin/micromamba shell hook -s bash)" \
    && micromamba create -y -n phc-meep -c conda-forge \
      python=3.11 \
      meep \
      numpy scipy matplotlib pandas \
      h5py mpi4py \
    && micromamba clean -a -y

RUN eval "$(/home/ubuntu/.local/bin/micromamba shell hook -s bash)" \
    && micromamba activate phc-meep \
    && python -c "import meep as mp; print('meep', getattr(mp, '__version__', 'unknown'))"

RUN mkdir -p /home/ubuntu/.cache/phc-meep \
    && echo ok > /home/ubuntu/.cache/phc-meep/.installed

RUN printf '%s\n' \
    'eval "$($HOME/.local/bin/micromamba shell hook -s bash)"' \
    'micromamba activate phc-meep' \
    >> /home/ubuntu/.bashrc
