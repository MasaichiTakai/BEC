#!/bin/bash

for v in `seq 0.03 0.03 1.500001`
do
    echo python run.py ${v} 1.5 ">>" logfile.txt >> joblist.sh
done

for v in `seq 0.002 0.002 0.2000001`
do
    echo python run.py ${v} 1.7 ">>" logfile.txt >> joblist.sh
done
