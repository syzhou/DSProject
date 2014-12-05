#!/bin/bash

./kmeans_2d_sequential input2d.csv output2dSeq.txt 2
mpirun -np 1 kmeans2DParallel input2d.csv output2dPara_1.txt 2
mpirun -np 2 kmeans2DParallel input2d.csv output2dPara_2.txt 2

diff output2dSeq.txt output2dPara_1.txt
diff output2dSeq.txt output2dPara_2.txt

