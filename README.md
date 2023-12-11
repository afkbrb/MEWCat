# A Fast Exact Solver with Theoretical Analysis for the Maximum Edge-Weighted Clique Problem

> Lu Liu, Mingyu Xiao, and Yi Zhou. A Fast Exact Solver with Theoretical Analysis for the Maximum Edge-Weighted Clique Problem. AAAI-24.

This repository contains the Appendix, code, and dataset for MEWCat.

The code of MEWCat is in mewcat.cpp.
Other .cpp and .h files are the code of the PLS heuristic algorithm implemented by Satoshi Shimizu. 

To compile, simply type "make" or use the following command:
g++ -O3 -o mewcat mewcat.cpp  weighted_graph.cpp phased_local_search.cpp phased_local_search_adjmatrix.cpp phased_local_search_adjlist.cpp

Then type "./mewcat <dataset>" to run MEWCat on the given dataset. For example, try this:
./mewcat dataset/dimacs/brock200_1
