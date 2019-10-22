#!/bin/bash

name=md

mkdir -p run01-trajectory

# -------------- clean up  --------------- #

rm -f run01-trajectory/*

# ---------------------------------------- #

params=aenet.prm
steps=1000000
timestep=2.0       # fs
writeout=0.500     # ps
MDtype=2           # NVT
temperature=800    # K

time ~/bin/dynamic.x-aenet-gfortran-mkl -k ${name}.key $name $params $steps $timestep \
               $writeout $MDtype $temperature > dynamic.out

# ---------------------------------------- #

mv *.out ${name}.*  run01-trajectory/
mv run01-trajectory/${name}.key .
mv run01-trajectory/${name}.xyz .
mv run01-trajectory/${name}.dyn .

# ---------------------------------------- #
