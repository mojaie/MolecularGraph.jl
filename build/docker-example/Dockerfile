FROM julia:1.11.5-bookworm AS jl
RUN apt-get update && apt-get install -y build-essential
WORKDIR /julia
COPY . .
RUN mkdir depot && mkdir /opt/julia
ENV JULIA_DEPOT_PATH=./depot
RUN julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); using MolecularGraph; write(stdout, normpath(joinpath(dirname(pathof(MolecularGraph)), "..")))' >> /etc/mgpath
RUN cd $(cat /etc/mgpath)"build" && make
CMD ["/bin/bash"]
