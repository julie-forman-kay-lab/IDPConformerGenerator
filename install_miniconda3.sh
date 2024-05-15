#!/bin/bash

sleep 1
mkdir -p miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_24.3.0-0-Linux-x86_64.sh -O miniconda3/miniconda.sh
sleep 1
bash miniconda3/miniconda.sh -b -u -p miniconda3
sleep 1
rm -rf miniconda3/miniconda.sh
