#!/bin/sh
#PBS -l walltime=24:00:00,nodes=1:ppn=4:gpus=1,vmem=32gb -j oe

source $HOME/init.sh
cd $PBS_O_WORKDIR

SEED=$PBS_ARRAYID

libbi sample @config.conf @posterior.conf