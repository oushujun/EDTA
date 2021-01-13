## adapt from https://github.com/Kapeel/edta_docker/blob/master/Dockerfile
FROM continuumio/miniconda3

RUN conda config --add channels defaults
RUN conda config --add channels bioconda
RUN conda config --add channels conda-forge
RUN conda config --add channels r

RUN conda create -n EDTA
RUN conda install -c bioconda -c conda-forge edta
RUN conda install tensorflow 'h5py<3.0.0'

RUN git clone https://github.com/oushujun/EDTA.git
RUN echo "source activate EDTA" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH

ENTRYPOINT [ "/EDTA/EDTA.pl" ]
