#!/bin/bash

conda env update --file requirements.yml --prune
sleep 1
python setup.py develop --no-deps
sleep 1
