mpi
===
This repository contains c files which use mpi(Message Passing Interface) to perform a variety of parallel computing operations:

forcecopy.c : Generates 'n' random particles and uses parallel computing to calculate 
              the coloumb-force on each particle due to all other particles. The
              process of computing the force scales with increased number of cores

matrixcopy.c: Generates 2 matrices of size 'n' x 'n' and calculates the product using 
              parallel computing. The process of computing the product scales with 
              increased number of cores
              
hamiltoncircle.c: For a given hamilton circle of 32 nodes with 16 connections between 
              the node apart from adjacent connections, parallel computing is used to 
              decide the configurations of the 16 connections to achieve the lowest 
              Average Path Lenght(A) and Diameter(D)
