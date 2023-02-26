# Inference-via-Message-Passing-IMP-for-Matrix-Completion
Tool for reproduce  Figure 3  Matrix Completion (Recomend System) using LDPC code Message Passing from article  B. -H. Kim, A. Yedla and H. D. Pfister, "IMP: A message-passing algorithm for matrix completion," 2010 6th International Symposium on Turbo Codes & Iterative Information Processing, Brest, France, 2010, pp. 462-466, doi: 10.1109/ISTC.2010.5613803.


https://arxiv.org/abs/1007.0481



in "IMP: A Message-Passing Algorithm for Matrix Completion"

1.  Extract IMP.zip 
2.  Open Matlab command window
3.  Type "mex generate_graph.c" and "mex mpa_bk.c" to compile mex files

A. To generate rmse curve in left panel (on Netflix Data Matrix 1)
4. Run "testing1" to test 
5. Run "plotting1" to generate a curve 

B. To generate rmse curve in right panel (on Netflix Data Matrix 2)
6. Run "testing2" to test 
7. Run "plotting2" to generate a curve 

Note: 
For us, it takes about 30~40 mins to generate each curve. To see rough shapes 
of resulting curves faster, changing a variable in plotting1.m/plotting2.m 
from "num_simulations = 20" to "num_simulations = 2" is suggested.



