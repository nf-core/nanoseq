FROM nfcore/base:dev
LABEL authors="Chelsea Sawyer" \
      description="Docker image containing all software requirements for the nf-core/nanoseq pipeline"

## THIS DOCKER FILE ISNT REQUIRED FOR PIPELINE
## ALL DOCKER CONTAINERS ARE CURRENTLY OBTAINED FROM EXTERNAL SOURCES
## WAITING FOR LINT TESTS TO BE UPDATED BEFORE DELETING IT

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-nanoseq-1.0dev/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-nanoseq-1.0dev > nf-core-nanoseq-1.0dev.yml
