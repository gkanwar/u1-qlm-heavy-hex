#!/bin/bash
# 16T29P10O final runs

shape=16T29P10O
i=7245385

geom_file=geoms/${shape}_geom.dat
init_file=geoms/${shape}_init.dat

for Kp_T in "0.40 0.73" "0.70 0.52" "2.00 0.55"; do
    Kp=$(echo ${Kp_T} | cut -f1 -d' ')
    T=$(echo ${Kp_T} | cut -f2 -d' ')
    echo Kp=${Kp} T=${T}
    for NT in 16 32 64 128; do
        dt=$(echo "(1/${T})*(1/${NT})" | bc -l)
        echo dt=${dt}
        ./bin/cluster.${NT}_8_4.release -t ${dt} -p ${Kp} -e 0.0 \
                      -f data_new/${shape}/final_KP${Kp}_T${NT} \
                      -s 100000 -m 100 -i 1000000  -r ${i} \
                      -b file -c file -y ${geom_file} -z ${init_file}
        ((i++))
    done
done
