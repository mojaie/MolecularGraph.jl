FROM julia:1.13.2-trixie AS builder
RUN apt-get update \
    && apt-get install -y --no-install-recommends build-essential make \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /julia
COPY . .
ENV JULIA_DEPOT_PATH=/julia/depot:
ENV PATH="/julia/depot/bin:${PATH}"
RUN mkdir -p /opt/julia
RUN julia -e 'using Pkg; Pkg.Apps.add("JuliaC"); Pkg.activate("."); Pkg.instantiate()'
RUN cd "$(julia --project=. -e 'using MolecularGraph; print(pkgdir(MolecularGraph))')" \
    && make build

FROM debian:trixie-slim
COPY --from=builder /opt/julia/libmolgraphjl/ /opt/julia/libmolgraphjl/
ENV LD_LIBRARY_PATH=/opt/julia/libmolgraphjl/lib
CMD ["/bin/sh"]
