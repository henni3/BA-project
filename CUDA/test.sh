#!/bin/bash
# this is inspired with help from chatGPT
# the files we need to run
paths=(swiss42.tsp berlin52.tsp kroA100.tsp kroA200.tsp kroB100.tsp kroC100.tsp)

for path in "${paths[@]}"
do

    ./test_gb ../Data/$path 1000
done

