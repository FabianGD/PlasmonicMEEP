#! /usr/bin/env bash

A=( $@ )

dir=${A[@]:0:1}
args=${A[@]:0}

cmd="ls $args"
echo $'Issuing the following cmd to SLURM:\n > '$cmd $'\n\n'
echo $($cmd)