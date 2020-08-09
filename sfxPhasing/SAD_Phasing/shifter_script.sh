#!/usr/bin/env bash
#SBATCH --image=docker:registry.services.nersc.gov/annagian/sfx_6
#SBATCH --nodes=1
#SBATCH --tasks-per-node=32
#SBATCH --constraint=haswell
#SBATCH --qos=debug
#SBATCH --gres=craynetwork:4
#SBATCH --job-name=sfx
#SBATCH --account=m3506
#SBATCH --mail-user=agiannakou@lbl.gov
#SBATCH --mail-type=ALL
#SBATCH --time=00:30:00
#run the application:
python batch_sub.py -rfl a2a.mtz -seq a2a.fasta -SFAC S -n 1 --shifter 1
