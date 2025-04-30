#!/usr/bin/bash

# usage:
# $ bash ipynb_to_html.sh 03-chrombpnet/01-train_models/03b-plot_correlation.ipyn

# convert notebook to HTML
# use the pretty-jupyter template https://github.com/JanPalasek/pretty-jupyter
jupyter nbconvert $1 --to html --template pj