FROM julia:1.12.1-bookworm

LABEL org.opencontainers.image.title="MOOSE"
LABEL org.opencontainers.image.description="Julia container for Mock Observation Of Synchrotron Emission"

ENV JULIA_DEPOT_PATH=/opt/julia-depot: \
    JULIA_PKG_PRECOMPILE_AUTO=0

WORKDIR /app

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        git \
        python3 \
    && rm -rf /var/lib/apt/lists/*

COPY Project.toml Manifest.toml ./

RUN julia --startup-file=no --project=/app -e 'using Pkg; Pkg.instantiate()'

COPY . .

RUN julia --startup-file=no --project=/app -e 'using Pkg; Pkg.precompile(); using MOOSE'

CMD ["julia", "--startup-file=no", "--project=/app", "-e", "using MOOSE; run_moose(help=true)"]
