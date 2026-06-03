#!/bin/bash
# Activation wrapper used as the container's ENTRYPOINT.
# Reads the env name from NXF_ENV (set by Nextflow via containerOptions)
# and activates the corresponding micromamba env before executing the command.

ENV_NAME="${NXF_ENV:-retro-base}"

eval "$(micromamba shell hook --shell bash)"
micromamba activate "$ENV_NAME"

exec "$@"
