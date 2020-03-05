FROM nfcore/base:1.9
LABEL authors="Chelsea Sawyer" \
      description="Docker image containing all software requirements for the nf-core/nanoseq pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a

# Add conda installation dir to PATH (instead of doing 'conda activate')
ENV PATH /opt/conda/envs/nf-core-nanoseq-1.0.0/bin:$PATH

# Dump the details of the installed packages to a file for posterity
RUN conda env export --name nf-core-nanoseq-1.0.0 > nf-core-nanoseq-1.0.0.yml
