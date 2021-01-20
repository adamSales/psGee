#!/bin/bash

#SBATCH -p short
#SBATCH -N 1
#SBATCH -n 50
#SBATCH -t 20:00:00
#SBATCH --mem 100G


R --no-save < simTot.r
