#!/bin/bash
b0=0.01  # beta = 1/kT
b=0.01
J=1
h=0
echo H C M X > ${1}.csv
for(( i=1;i<=100;i++ ))
do
    b=`echo "${b0} * $i"| bc -l`
    ./ising.run ${b} ${J} ${h} >> ${1}.csv
    echo "$i / 100"
done
