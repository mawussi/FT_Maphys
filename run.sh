#!/bin/bash

make -j 20
make -j 20 examplelv
cd examples
# run young1c with Fgmres
mpirun -np 4 ./zmph_examplelv complex_young1c.lvin 
#create tags 
find . -name "*.F90"   -print > taglist.txt
find . -name "*.[Ff]"   -print >> taglist.txt
cat taglist.txt | etags -
#mpirun -np 4 ./dmph_examplelv real_bcsstk17.lvin > dmph_bcsstk17.lvout #&& more  dmph_bcsstk17.lvout