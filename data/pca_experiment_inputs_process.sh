#!/bin/bash
# ls experiment_inputs | grep 'pos' | cut -c-10 > experiment_inputs/go_list.txt
for f in pca_experiment_inputs/*.txt
do
    awk '{printf("%010d %s\n", NR, $0)}' ${f} > ${f}.txt 
done

