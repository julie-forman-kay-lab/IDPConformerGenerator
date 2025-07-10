FROM ubuntu:22.04

RUN apt update
RUN apt install python3 python3-pip -y

# This part is copied from DSSP
# https://github.com/cmbi/dssp/blob/master/Dockerfile
RUN apt-get install -y make rsync wget
RUN apt-get install -y git g++ libboost-all-dev libbz2-dev doxygen xsltproc docbook docbook-xsl docbook-xml autoconf automake autotools-dev

RUN mkdir -p /deps

# Install libzeep
RUN git clone https://github.com/mhekkel/libzeep.git /deps/libzeep ;\
    cd /deps/libzeep ;\
    git checkout tags/v3.0.3
# XXX: Workaround due to bug in libzeep's makefile
RUN sed -i '71s/.*/\t\$\(CXX\) \-shared \-o \$@ \-Wl,\-soname=\$\(SO_NAME\) \$\(OBJECTS\) \$\(LDFLAGS\)/' /deps/libzeep/makefile
WORKDIR /deps/libzeep
# XXX: Run ldconfig manually to work around a bug in libzeep's makefile
RUN make ; make install ; ldconfig

WORKDIR /deps
RUN git clone https://github.com/cmbi/dssp dssp
WORKDIR /deps/dssp
RUN ./autogen.sh && ./configure && make && make install

RUN mkdir -p /app
WORKDIR /app

RUN git clone https://github.com/julie-forman-kay-lab/IDPConformerGenerator idpconfgen
WORKDIR /app/idpconfgen

RUN ./install_miniconda3.sh
RUN ./miniconda3/bin/pip install -r requirements.txt
RUN ./miniconda3/bin/python setup.py develop --no-deps

RUN echo 'export PATH="/app/idpconfgen/miniconda3/bin/:$PATH"' >> ~/.bashrc

WORKDIR /deps
RUN git clone https://github.com/THGLab/MCSCE
WORKDIR /deps/MCSCE
RUN /app/idpconfgen/miniconda3/bin/pip install tensorflow tdqm pathos
RUN /app/idpconfgen/miniconda3/bin/python setup.py develop --no-deps

WORKDIR /app/idpconfgen
RUN ./miniconda3/bin/python setup.py develop --no-deps

