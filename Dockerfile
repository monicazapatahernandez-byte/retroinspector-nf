FROM condaforge/miniforge3:24.7.1-0

WORKDIR /opt/retroinspector-nf

COPY env.yaml /tmp/env.yaml
COPY env_sniffles2.yaml /tmp/env_sniffles2.yaml
COPY env_repeatmasker_t2t.yaml /tmp/env_repeatmasker_t2t.yaml
COPY r.yaml /tmp/r.yaml

RUN mamba env create -n retro-base -f /tmp/env.yaml && \
    mamba env create -n retro-sniffles2 -f /tmp/env_sniffles2.yaml && \
    mamba env create -n retro-rm-t2t -f /tmp/env_repeatmasker_t2t.yaml && \
    mamba env create -n retro-r -f /tmp/r.yaml && \
    mamba clean -a -y

ENV PATH=/opt/conda/envs/retro-base/bin:$PATH
ENV RETRO_BASE_ENV=/opt/conda/envs/retro-base
ENV RETRO_SNIFFLES2_ENV=/opt/conda/envs/retro-sniffles2
ENV RETRO_RM_T2T_ENV=/opt/conda/envs/retro-rm-t2t
ENV RETRO_R_ENV=/opt/conda/envs/retro-r
