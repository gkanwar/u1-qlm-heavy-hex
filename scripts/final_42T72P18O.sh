#!/bin/bash
# 42T72P18O final runs

shape=42T72P18O
i=811167385

geom_file=geoms/${shape}_geom.dat
init_file=geoms/${shape}_init.dat

for Kp_T in "0.40 0.72" "0.70 0.48" "2.00 0.59"; do
    Kp=$(echo ${Kp_T} | cut -f1 -d' ')
    T=$(echo ${Kp_T} | cut -f2 -d' ')
    echo Kp=${Kp} T=${T}
    for NT in 16 32 64 128; do
        dt=$(echo "(1/${T})*(1/${NT})" | bc -l)
        echo dt=${dt}
        ./bin/cluster.${NT}_8_5.release -t ${dt} -p ${Kp} -e 0.0 \
                      -f data_new/${shape}/final_KP${Kp}_T${NT} \
                      -s 100000 -m 100 -i 1000000  -r ${i} \
                      -b file -c file -y ${geom_file} -z ${init_file}
        ((i++))
    done
done
