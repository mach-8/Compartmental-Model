#!/bin/bash

# base directory
base_dir="$( cd "$( dirname "$0" )" && pwd )"

# output directory
output_dir="$base_dir/OUTPUTS"

# input directory
input_dir="$base_dir/INPUTS"

# input file
input_file="$input_dir/input_parameters.yaml"

# log file
log_file="$output_dir/simulator_log.txt"

# go into the inputs directory and generate input file
cd INPUTS/
./input_parameters.py
cd ..


# clean up output directory
mkdir OUTPUTS/
cd OUTPUTS/
rm *.dat *.txt
cd ..

# clean up figure directory
cd FIGURES/
rm *.pdf
cd ..

# compile source files
cd SOURCES/
make clean
make all

# run executable
./simulator $input_file $output_dir >> $log_file
cd ..

# create figures
cd FIGURES/
./createGraphs.plt
cd ..