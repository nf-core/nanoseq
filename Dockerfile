FROM nfcore/base:1.7
LABEL authors="Chelsea Sawyer" \
      description="Docker image containing all requirements for nf-core/nanodemux pipeline"

COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
ENV PATH /opt/conda/envs/nf-core-nanodemux-1.0dev/bin:$PATH
