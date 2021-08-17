FROM ubuntu:20.04
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update

COPY . /blimpy
WORKDIR /blimpy

RUN cat dependencies.txt | xargs -n 1 apt install --no-install-recommends -y

RUN cd tests && bash download_data.sh && cd ..

RUN python3 -m pip install --user -r requirements.txt
RUN python3 setup.py install --user
RUN python3 -m pip install --user -r requirements_test.txt
RUN python3 setup.py test

RUN rm -fr tests/test_data
RUN find . -path '*/__pycache__*' -delete

WORKDIR /home
