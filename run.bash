#!/usr/bin/env bash

for g in 0.0035 0.0105 0.021 0.035 0.105
do
    for syn in sphere random5
    do
               python main.py    -g=$g -s=$syn
               python main.py -c -g=$g -s=$syn
    done
done
