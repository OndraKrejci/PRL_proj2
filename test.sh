#!/bin/bash

if [ $# -lt 1 ]; then
	echo "Requires one argument representing nodes of a binary tree as an array"
	exit 1
fi

n=${#1};
let procs=(2*n)-2;
if [ $procs -eq 0 ]; then
	procs=1;
fi

mpic++ --prefix /usr/local/share/OpenMPI -o pro pro.cpp

mpirun --prefix /usr/local/share/OpenMPI --oversubscribe -np $procs pro $1

rm -f pro
