#!/bin/sh
#PBS -l walltime=48:00:00,nodes=1:ppn=16:gpus=1,vmem=64gb -j oe

source $HOME/init.sh
cd $PBS_O_WORKDIR

libbi sample @config.conf @posterior.conf --filter bootstrap --output-file results/posterior_bootstrap.nc
