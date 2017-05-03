#!/bin/bash


#PBS -mabe
#PBS -l walltime=66:55:44
#PBS -N HEdose_fast_sim
#PBS -t 1-25


module load apps/plink-1.90
module load languages/R-3.2.0-ATLAS 

cd $PBS_O_WORKDIR

#nobs=6500
#nsnps=24000
nrep=12
#hsqA=0.01
#hsqP=0.2
#delta_Maf=0
#rel_Cor=0
npairvar=1
avg_LD=0

for delta_Maf in 0 0.05
do
    for rel_Cor in 0 0.25
    do
        for nsnps in 500.0 1000.0 3000.0
        do
            for nobs in 1000 2000 4000
            do
                hsqA=$(echo "scale=8;0.1 * ($nsnps / 3000.0)"| bc)
                hsqD=$(echo "scale=8;0.1 * ($nsnps / 3000.0)"| bc)
                hsqP=$(echo "scale=8;0.1 * ($nsnps / 3000.0)"| bc)
                Rscript --vanilla lrt_sim_HEDOSE_0.6.R ${PBS_ARRAYID} ${hsqA} ${hsqD} ${hsqP} ${delta_Maf} ${rel_Cor} ${avg_LD} ${npairvar} ${nobs} ${nsnps} ${nrep}
                
            done
        done
    done
done
