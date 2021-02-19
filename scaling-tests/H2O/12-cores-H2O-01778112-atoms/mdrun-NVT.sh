#!/bin/bash

name=md

mkdir -p run01-trajectory

# -------------- clean up  --------------- #

rm -f run01-trajectory/*

# ---------------------------------------- #

params=aenet.prm
steps=1
timestep=1.00         # fs
writeout=1.000       # ps
MDtype=2              # NVT
temperature=300       # K

dynamic_x="/moto/urban/users/na2782/Programs/tinker-8.7.1-aenet-v2.1.0/source/dynamic.x"

time ${dynamic_x} -k ${name}.key $name $params $steps $timestep \
          $writeout $MDtype $temperature > dynamic.out

# ---------------------------------------- #

mv *.out ${name}.*  run01-trajectory/
mv run01-trajectory/${name}.key .
mv run01-trajectory/${name}.xyz .
mv run01-trajectory/${name}.dyn .

# ---------------------------------------- #
