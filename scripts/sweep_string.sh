#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --output=slurm_logs/slurm-%j.out

### Sweep over KP for given string geometry

WORKDIR=/space2/kanwar/u1-qlm-heavy-hex
cd ${WORKDIR}

function usage() {
    echo "Usage: $0 <KP> <out_prefix> -- <extra args>"
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

KP=$1
out_prefix=$2
NROWS=${NROWS:-64}
NCOLS=${NCOLS:-64}
NT=${NT:-32}
srun ./bin/cluster.${NT}_${NROWS}_${NCOLS}.release -p ${KP} -f "${out_prefix}.KP${KP}" "${@:4}"

