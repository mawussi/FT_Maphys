# FT_MaPHyS

## FT_MaPHyS a fault tolerant MaPHyS
[**MaPHyS**](https://gitlab.inria.fr/solverstack/maphys/maphys/) is parallel linear solver couples direct and iterative
  approaches. The underlying idea is to apply to general unstructured linear systems domain decomposition
  ideas developed for the solution of linear systems arising from PDEs. The interface problem,
  associated with the so called Schur complement system, is solved using a block preconditioner
  with overlap between the blocks that is referred to as Algebraic Additive Schwarz. To better exploit the architecture of
  emerging parallel platforms the solver exploits two levels of parallelism
  (between the blocks and within the treatment of the blocks). This enables us to exploit a large number
  of cores with a moderate number of blocks which ensures a reasonable convergence behaviour especially for non SPD linear systems.For these latter problems a two-level preconditionner, based on Geneo ideas for the coarse grid construction, is available through various implementations.
  
## Fault simulation in FT_MaPHyS
 We simulate a process fault by overwriting its
dynamic data (an actual fault injection is an orthogonal problem out
of the present scope) and we then use either data redundancy or [numerical 
interpolation](https://hal.inria.fr/hal-01323192/file/final_nlaa.pdf) 
techniques (or both) to regenerate the lost dynamic data. We simulate
the crash of one single process (denoted *single process fault*) and 
the crash of multiple concurrent processes that are
neighbors with respect to the domain decomposition (denoted
*neighbor processes fault case*). When a single fault occurs, we
exploit data redundancy intrinsic to **MaPHyS** to retrieve all lost
dynamic data. When faults are simultaneously injected into neighbor
processes, part of the data is definitely lost on all processes; the
strategy then consists in exploiting data redundancy wherever
possible enhanced with an numerical interpolation scheme to regenerate definitely lost
dynamic data.

## Requirements

#### Graph partitioning 
* **Metis**: available at [metis-5.1.0.tar.gz](http://glaros.dtc.umn.edu/gkhome/fetch/sw/metis/metis-5.1.0.tar.gz).
* **Scotch**: available at [scotch_6.0.4.tar.gz](http://gforge.inria.fr/frs/download.php/file/34618/scotch_6.0.4.tar.gz).

#### Sparse direct solver
* **Mumps** for sparse LU factorization available at [MUMPS_5.0.1.tar.gz](http://mumps.enseeiht.fr/MUMPS_5.0.1.tar.gz)
* **PaSTiX** an alternative for Mumps available at   [pastix-6.0.0.tar.gz](http://pastix.gforge.inria.fr/files/README-txt.html)

#### BLAS & LAPACK
* **MKL** is now free for academics (students and researchers) available at [https://software.intel.com/en-us/articles/free-mkl](https://software.intel.com/en-us/articles/free-mkl)

OR
* **Netlib BLAS** no optimized BLAS routines, available at [BLAS-3.6.0.tgz](http://www.netlib.org/blas/blas-3.6.0.tgz)
* **Netlib LAPACK** no optimized LAPACK routines, available at [LAPACK-3.6.1.tgz](http://www.netlib.org/lapack/lapack-3.6.1.tgz) 
 
#### Hardware locality
* **Hwloc** for more information on the hardware topology,  available at [hwloc-1.11.3.tar.gz](https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.3.tar.gz)

## Installation 
1. git clone https://github.com/mawussi/FT_Maphys.git
2. cd FT_Maphys
3. configuration in make.inc, it consists of providing path to the required libraries 
4. make and make install
