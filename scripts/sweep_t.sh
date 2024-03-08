#!/bin/bash

### Sweep over NT and dt to allow zero-temperature, time-continuum extrapolation

function usage() {
    echo "Usage: $0 <base_dt> <out_prefix> -- <extra args>"
}

if [[ "$1" == "" ]]; then
    usage
    exit 1
fi
if [[ "$2" == "" ]]; then
    usage
    exit 1
fi
if [[ "$3" != "--" ]]; then
    usage
    exit 1
fi

base_dt=$1
out_prefix=$2
NROWS=2
NCOLS=2
for NT in 64 80 96 128; do
    dt=$(echo "scale=4; ${base_dt}*${NT}/64" | bc)
    ./bin/cluster.${NT}_${NROWS}_${NCOLS}.release \
                  -t "${dt}" -f "${out_prefix}.T${NT}_dt${dt}" \
                  "${@:4}" || exit 1
done
