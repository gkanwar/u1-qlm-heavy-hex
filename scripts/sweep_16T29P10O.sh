#!/bin/bash
# 16T29P10O lattice sweep

shape=16T29P10O
i=7876432
geom_file=geoms/${shape}_geom.dat
init_file=geoms/${shape}_init.dat
for dt in $(seq 0.01 0.01 0.05); do
    for T in 16 32 64 128; do
        for Kp in 0.40 0.70 2.00; do
            ./bin/cluster.${T}_8_4.release -t ${dt} -p ${Kp} -e 0.0 \
                          -f data_new/${shape}/KP1.00_T${T}_dt${dt} \
                          -s 100000 -m 100 -i 1000000  -r ${i} \
                          -b file -c file -y ${geom_file} -z ${init_file}
            ((i++))
        done
    done
done
