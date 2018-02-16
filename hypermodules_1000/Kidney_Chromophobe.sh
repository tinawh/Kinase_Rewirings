#!/bin/bash 
module load java/1.6.0_21 
        java -jar /.mounts/labs/reimandlab/private/users/thuang/HyperModules_1.0.2_CMD/HyperModulesCMD-1.0.2.jar -n /.mounts/labs/reimandlab/private/users/thuang/data_2/read-only/network_interaction_data_hugo.tsv -s /.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/mutation_data/mut_Kidney_Chromophobe_hugo.csv -c /.mounts/labs/reimandlab/private/users/thuang/data_2/01-08-18/clinical_data/Kidney_Chromophobe_clin.csv -S 1000 -t logrank -p 0.05 -C 32 > /.mounts/labs/reimandlab/private/users/thuang/data_2/hypermodules_1000/hypermodules_1000_Kidney_Chromophobe.txt