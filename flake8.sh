#!/bin/bash

# W293 blank line contains whitespaces
# W503 line break before binary operator, actually favoured by PEP8
# -hang-closing, allows:
#my_func(
#    var1,
#    var2,
#    )
#E402,F401,W605
clear
flake8 --hang-closing --ignore=W293,W503 $1
