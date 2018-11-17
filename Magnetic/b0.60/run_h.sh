#!/bin/bash
h0=0.005  # beta = 1/kT
b=${2}
J=1
h=${h0}
echo H C M X > ${1}${2}.csv
./ising.run ${b} ${J} 0 >> ${1}${2}.csv
for(( i=1;i<=50;i++ ))
do
    ./ising.run ${b} ${J} ${h} >> ${1}${2}.csv
    h=`echo "${h0} * $i"| bc -l`
    echo "$i / 50"
done

python3 ../../plt2.py ${1}${2}.csv ${h0}
