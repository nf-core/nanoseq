FROM nfcore/base:1.7
LABEL authors="Chelsea Sawyer" \
      description="Docker image containing all requirements for nf-core/nanoseq pipeline"

## THIS DOCKER FILE ISNT REQUIRED FOR PIPELINE
## ALL DOCKER CONTAINERS ARE CURRENTLY OBTAINED FROM EXTERNAL SOURCES
## WAITING FOR LINT TESTS TO BE UPDATED BEFORE DELETING IT

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-nanoseq-1.0dev > nf-core-nanoseq-1.0dev.yml
ENV PATH /opt/conda/envs/nf-core-nanoseq-1.0dev/bin:$PATH
