#!/usr/bin/env bash

cd /.mounts/labs/reimandlab/private/users/thuang/bin_2/hypermodules_100/

for file in *.sh; do 
	qsub -cwd -V -l h_vmem=70G $file
done