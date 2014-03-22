mpi
===
This repository contains c files which use mpi to perform a variety of operations:

forcecopy.c : Generates 'n' random particles and uses parallel computing to calculate 
              the coloumb-force on each particle due to all other particles

matrixcopy.c: Generates 2 matrices of size 'n' x 'n' and calculates the product using 
              parallel computing
