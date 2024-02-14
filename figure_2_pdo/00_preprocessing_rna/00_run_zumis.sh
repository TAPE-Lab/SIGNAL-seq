#!/bin/bash -l

# Run the zUMIs pipeline
# Each sublibrary was run via an independant .yaml file
# -c specifies the use of a conda environment
~/zUMIs/zUMIs.sh -c -y ~/PATH/ex0015_rna_76.yaml
