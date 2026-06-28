ARG JULIA_VERSION=1.12.6
FROM julia:${JULIA_VERSION}-bookworm

LABEL org.opencontainers.image.title="MOOSE"
LABEL org.opencontainers.image.description="Julia container for Mock Observation Of Synchrotron Emission"

ENV JULIA_DEPOT_PATH=/opt/julia-depot: \
    JULIA_PKG_PRECOMPILE_AUTO=0 \
    JULIA_NUM_THREADS=auto \
    PYTHONUNBUFFERED=1 \
    MOOSE_PROJECT=/app

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        fontconfig \
        git \
        libgl1 \
        libx11-6 \
        libxext6 \
        libxrender1 \
        libxt6 \
        python3 \
        tini \
    && rm -rf /var/lib/apt/lists/*

COPY Project.toml Manifest.toml ./

RUN julia --startup-file=no --project=/app -e 'using Pkg; Pkg.instantiate()'

COPY README.md Version.toml setup.jl ./
COPY config ./config
COPY python ./python
COPY src ./src
COPY test ./test

RUN julia --startup-file=no --project=/app -e 'using Pkg; Pkg.precompile(); using MOOSE'

RUN mkdir -p /data

VOLUME ["/data"]

ENTRYPOINT ["tini", "--"]

CMD ["julia", "--startup-file=no", "--project=/app", "-e", "using MOOSE; run_moose(help=true)"]
