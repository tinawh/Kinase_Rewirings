#!/usr/bin/python3
import os

def write_hypermodules_script(mutation_dir, clinical_dir, file_location = "/.mounts/labs/reimandlab/private/users/thuang/bin_2/hypermodules_biogrid_100/", result_location = "/.mounts/labs/reimandlab/private/users/thuang/data_2/hypermodules_biogrid_100/"):
        if not os.path.exists(file_location):
                os.makedirs(file_location)

        if not os.path.exists(result_location):
                os.makedirs(result_location)

        mut_files = sorted(os.listdir(mutation_dir))
        clin_files = sorted(os.listdir(clinical_dir))

        muts = [x[4:-9] for x in mut_files]
        clins = [x[:-9] for x in clin_files]
        assert muts == clins, "mutation cancer types different from clinical"

        lst = []
        for a in mut_files:
                lst.append(a[4:-9])

        file_names = []
        for name in lst:
                file_names.append(file_location + name + ".sh")

        files  = ""
        for i in range(len(mut_files)):
                files += """java -jar /.mounts/labs/reimandlab/private/users/thuang/HyperModules_1.0.2_CMD/HyperModulesCMD-1.0.2.jar -n /.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/biogrid_network.tsv -s {0}{1} -c {2}{3} -S 100 -t logrank -p 0.05 -C 32 > {4}hypermodules_100_{5}.txt
        """.format(mutation_dir, mut_files[i], clinical_dir, clin_files[i], result_location, lst[i])
        
        commands = files.split("\n")
        commands = commands[:-1]

        for i in range(len(commands)):
                with open(file_names[i], "w+") as file:
                        file.write("#!/bin/bash \nmodule load java/1.6.0_21 \n")
                        file.write(commands[i])
        print("done")

if __name__ == "__main__":
        mutation_dir = "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/mutation_data/"
        clinical_dir = "/.mounts/labs/reimandlab/private/users/thuang/data_2/02-19-18/clinical_data/"
        # network_dir = '/Users/t2huang/Desktop/network_interaction_data_hugo.tsv'
        write_hypermodules_script(mutation_dir, clinical_dir)
