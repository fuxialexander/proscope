FROM mambaorg/micromamba

WORKDIR /app

RUN micromamba install -n base -c conda-forge -c bioconda -y python=3.10 git pip biopython nglview tqdm matplotlib pandas xmlschema seaborn numpy py3Dmol 
ARG MAMBA_DOCKERFILE_ACTIVATE=1


USER root
RUN mkdir /data
RUN apt-get update && apt-get install -y git ssh && apt-get clean && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

COPY --chown=$MAMBA_USER:$MAMBA_USER modules /app/modules
COPY --chown=$MAMBA_USER:$MAMBA_USER app /app/app

RUN cd modules/proscope && pip3 install .

WORKDIR /app


