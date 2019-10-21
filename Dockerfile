FROM nfcore/base:dev
LABEL authors="Chelsea Sawyer" \
      description="Docker image containing all requirements for nf-core/nanoseq pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
RUN conda env export --name nf-core-nanoseq-1.0dev > nf-core-nanoseq-1.0dev.yml
ENV PATH /opt/conda/envs/nf-core-nanoseq-1.0dev/bin:$PATH
