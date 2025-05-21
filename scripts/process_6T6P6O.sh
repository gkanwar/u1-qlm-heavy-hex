#!/bin/bash

for dt in $(seq 0.01 0.01 0.05); do
    for T in 16 32 64 128; do
        python analyze.py --prefix data_new/6T6P6O/KP1.00_T${T}_dt${dt}
    done
done
