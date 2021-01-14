FROM continuumio/miniconda3

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge

RUN conda install -c conda-forge mamba && conda clean --all -f --yes
RUN mamba install -c conda-forge -c bioconda edta python=3.6 tensorflow=1.14 'h5py<3' && mamba clean --all -f --yes

RUN git clone https://github.com/oushujun/EDTA.git
ENV PATH /opt/conda/envs/env/bin:$PATH

