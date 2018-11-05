# Base Image
FROM python:2.7.15-stretch

# Metadata
LABEL base.image="pipeline-tools:latest"
LABEL version="1"
LABEL software="Pipeline-Tools"
LABEL software.version="latest"
LABEL description="Python package for various NGS utility functions."
LABEL tags="recode vcf genotype"

# Maintainer
MAINTAINER Dave Lab <lab.dave@gmail.com>

# update the OS related packages
RUN apt-get update

# install required dependencies for QCParser
RUN pip install numpy
RUN pip install scipy
RUN pip install pandas
RUN pip install matplotlib
RUN pip install pyVCF

# make directory to store tools such as QCParser
RUN mkdir tools

# get the QCParser from GitHub
RUN git clone --branch master https://github.com/labdave/Pipeline-Tools.git /tools/PipelineTools

RUN chmod 777 -R /tools/PipelineTools

ENV PATH /tools/PipelineTools:$PATH

CMD ["python", "RecodeVCF.py"]
