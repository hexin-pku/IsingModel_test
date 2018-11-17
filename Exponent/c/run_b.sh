#!/bin/bash

# this bash run a series Ising simulation in different tenperature
# it covers b=0~1, where phase transition occurs around about b=0.44

J=1
h=0
echo H C M X > ${1}.csv
for b in 0.42 0.422 0.424 0.426 0.428 \
0.43 0.432 0.434 0.436 0.438 \
0.44 0.442 0.444 0.446 0.448 \
0.45 0.452 0.454 0.456 0.458 0.46
do
    #b=`echo "${b0} * $i"| bc -l`
    ./ising.run ${b} ${J} ${h} >> ${1}.csv
    echo ${b}
done

