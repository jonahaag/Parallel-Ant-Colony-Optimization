# Parallel-Ant-Colony-Optimization

This repository contains the project work for the course "Parallel Computations for Large-Scale Problems" at the Royal Institute of Technology (KTH) in the spring term 2023, a joint work with Nolwenn Deschand.
The purpose of the project was to implement Parallel Ant Colony Optimization (Parallel ACO) to solve the Traveling Salesman Problem (TSP) using C and the Message Passing Interface (MPI).

ACO is a biologically inspired optimization metaheuristic, described in great detail in ["Ant Colony Optimization"](https://ieeexplore.ieee.org/abstract/document/4129846?casa_token=eaotD19GV74AAAAA:fQ9trQuIiu6cP0CncXA_piODpGBRPBkabmRdQQ4b68YK_oW-PvvXnPKIayo4rzgOKpBtyReMtxk) by Dorigo et al.
This implementation is based on ["Max-Min Ant System"](https://www.sciencedirect.com/science/article/pii/S0167739X00000431?casa_token=txj73Si8h2AAAAAA:n3MdoZZ1KiMjmPDpJEAS8N-7gxHv8e-9wJ0vLjWMk0jVcy3PI6F1xioLOQsEYH2ljIF0UGqJv_I) by Thomas Stützle and Holger Hoos and was tested using various examples from the [TSPLIB](http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/) of the University of Heidelberg.

The approach and results are described and analyzed in the [report](report_parallel_aco.pdf).
The experiments mentioned in the paper were performed on KTH's supercomputer Dardel at [PDC](https://www.pdc.kth.se).

To test, both the sequential version [`aco.c`](src/aco.c) and the parallel implementation [`aco_mpi.c`](src/aco_mpi.c) are provided in this repository with three example TSP instances of different size.

Exemplary workflow for `berlin52`:

```bash
cc -O3 -o src/aco src/aco.c
./aco data/berlin52/berlin52.csv data/berlin52/params_berlin52.txt
```

```bash
mpicc -O3 -o src/aco_mpi src/aco_mpi.c
mpirun -n 4 ./aco_mpi data/berlin52/berlin52.csv data/berlin52/params_berlin52.txt
```