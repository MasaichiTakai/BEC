#!/bin/bash

for v in `seq 0.002 0.002 0.100001`
do
    echo mkdir /Users/takaimasaichi/Desktop/a=1.6/v=${v} >> transportBEC.sh
    echo scp b37107@laurel.kudpc.kyoto-u.ac.jp:supercomBEC/v=${v}a=1.6/* /Users/takaimasaichi/Desktop/a=1.6/v=${v} >> transportBEC.sh
    
done
