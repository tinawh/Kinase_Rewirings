#!/usr/bin/env bash

cd /.mounts/labs/reimandlab/private/users/thuang/bin_2/hypermodules_biogrid_10/

# for file in *.sh; do 
qsub -cwd -V -l h_vmem=80G Kidney_renal_clear_cell_carcinoma.sh
qsub -cwd -V -l h_vmem=80G Lymphoid_Neoplasm_Diffuse_Large_B-cell_Lymphoma.sh
qsub -cwd -V -l h_vmem=80G Pancreatic_adenocarcinoma.sh

# done