#!/bin/bash
ls experiment_inputs_subset | grep 'pos' | cut -c-10 > experiment_input_go_list.txt
