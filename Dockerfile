FROM python:3.12.16-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/colinvwood/fq-manifestor"

RUN apt-get update \
    && apt-get install -y --no-install-recommends bash procps \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/fq-manifestor
COPY pyproject.toml README.md LICENSE ./
COPY fq_manifestor ./fq_manifestor
RUN python -m pip install --no-cache-dir typer==0.27.3 .
