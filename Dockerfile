
FROM mambaorg/micromamba:1.5.8

USER root

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps wget curl ca-certificates bash git unzip \
    && rm -rf /var/lib/apt/lists/*

COPY env.yaml                   /opt/envs/env.yaml
COPY env_sniffles2.yaml         /opt/envs/env_sniffles2.yaml
COPY r.yaml                     /opt/envs/r.yaml
COPY env_repeatmasker_t2t.yaml  /opt/envs/env_repeatmasker_t2t.yaml

RUN micromamba env create -y -n retro-base       -f /opt/envs/env.yaml \
 && micromamba env create -y -n retro-sniffles2  -f /opt/envs/env_sniffles2.yaml \
 && micromamba env create -y -n retro-rm-t2t     -f /opt/envs/env_repeatmasker_t2t.yaml \
 && micromamba env create -y -n retro-r          -f /opt/envs/r.yaml \
 && micromamba clean --all --yes

ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV PATH=/opt/conda/condabin:/opt/conda/bin:$PATH
 

RUN chmod -R a+rx /opt/conda

COPY activate-and-run.sh /usr/local/bin/activate-and-run.sh
RUN chmod +x /usr/local/bin/activate-and-run.sh
ENTRYPOINT ["/usr/local/bin/activate-and-run.sh"]

WORKDIR /work
