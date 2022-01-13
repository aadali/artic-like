#!/bin/bash

wd=`pwd`

echo export PATH=\$PATH:$wd/texlive >> ~/.bashrc    # source ~/.bashrc

condaDir=$(dirname $(dirname $(which conda)))

sed -i "/^prefix/c prefix: $condaDir/envs/artic-like" env/artic-like.env.yaml

echo "conda env create -f env/artic-like.env.yaml"

conda env create -f env/artic-like.env.yaml

##
