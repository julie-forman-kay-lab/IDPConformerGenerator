#!/bin/bash

conda env update --file requirements.yml
sleep 1
python setup.py develop --no-deps
sleep 1
