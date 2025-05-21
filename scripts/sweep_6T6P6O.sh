#!/bin/bash
# 6T6P6O lattice sweep

i=12389
for dt in $(seq 0.01 0.01 0.05); do
    for T in 16 32 64 128; do
        ./bin/cluster.${T}_2_2.release -t ${dt} -p 1.0 -e 0.0 -f data_new/6T6P6O/KP1.00_T${T}_dt${dt} -s 100000 -m 100 -i 1000000 -b obc -c cold -r ${i}
        ((i++))
    done
done
