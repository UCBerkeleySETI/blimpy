FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update
RUN apt-get install --no-install-recommends -y \
    python3-pip \
    python3-dev \
    gfortran \
    bitshuffle \
    curl \
    git

COPY . /blimpy
WORKDIR /blimpy

#RUN cd tests && bash download_data.sh && cd ..

RUN python3 -m pip install -r requirements.txt
RUN python3 -m pip install -e .[full]
RUN python3 setup.py install
RUN python3 setup.py test

WORKDIR /home
