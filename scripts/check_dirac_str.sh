#!/bin/bash

for r in 4 6; do
    ./bin/cluster.16_8_8.release \
        -t 0.05 -p 0.70 -e 0.0 \
        -f data_new/test_dirac_str/test_v1_r${r} \
        -i 1000 -s 1000 -m 1 -b pbc -c cold \
        -d geoms/8_8_r${r}_dirac_str.dat -r 1234 >/dev/null
    ./bin/cluster.16_8_8.release \
        -t 0.05 -p 0.70 -e 0.0 \
        -f data_new/test_dirac_str/test_v2_r${r} \
        -i 1000 -s 1000 -m 1 -b pbc -c file \
        -z geoms/8_8_r${r}_init2.dat \
        -d geoms/8_8_r${r}_dirac_str2.dat -r 1234 >/dev/null
    echo "Hx:"
    sha256sum data_new/test_dirac_str/test_v1_r${r}.Hx.dat
    sha256sum data_new/test_dirac_str/test_v2_r${r}.Hx.dat
    echo "Ex:"
    sha256sum data_new/test_dirac_str/test_v1_r${r}.Ex.dat
    sha256sum data_new/test_dirac_str/test_v2_r${r}.Ex.dat
done

