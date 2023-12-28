# Bachelor Project DIKU 2023 - Traveling Salesmens Problem (TSP)
Supervisor: Associate Professor, Cosmin Eugen Oance  
Students: Henriette Naledi Winther Hansen and Jóhann Utne

# Structure of Github Repository
- The _report.pdf_ document is our BSc thesis: it presents the description of the project including a literature survey on TSP-related work. In particular it presents several algorithms aimed at solving TSP, then it chooses from these the one most suitable for execution on GPU hardware, then it zooms on the parallelization (and optimization) strategy and finally it presents an experiment that systematically evaluates the accuracy and the runtime (GPU) performance of our implementations.
- The _CUDA directory_ contains the CUDA implementation of the 2-Opt algorithm, together with publicly available code from O’Neil et.al (https://userweb.cs.txstate.edu/~mb92/papers/pdpta11b.pdf).   
- The _Data directory_, contains all the datasets acquired from the TSBLIB95 library (http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/).  
- The _Futhark directory_ contains the Futhark implementation of the 2-Opt algorithm.  
- The _Sequentielt directory_, contains the sequential implementations of the 2-Opt, ACO, and SA algorithms.  

# How to Run the Project
To compile and run the _CUDA implementation_ with a sample program go to the CUDA folder and insert,  
```
make prepro && ./prepro 
```
to terminal. This program will run on the dataset KroA100.tsp with 75000 climbers.  

To compile the _Futhark implementation_ go to the Futhark folder and insert
```
futhark cuda TSP.fut
```
to terminal. To run insert
```
./TSP −t /dev/stderr −r 10 < ../Data/berlin52.txt 
```
this will print the time 10 times and the final number (ending with i32) will be the best cost of the
tour computed. When loading the dataset berlin52.txt it will automatically run with 100 climbers.  

To compile the _sequential implementation_ go to the Seqventielt folder and insert
```
make all
```
in terminal. To run the programs different arguments are needed. To run the _ACO_ program insert
```
./aco <../Data/dataset.tsp> <number of ants> < number of iterations>
```
in terminal with a dataset, a given number of ants and a given number of iterations. To run the _SA_ program insert
```
./sa <../Data/dataset.tsp>
```
in terminal with a dataset. To run the _2-Opt_ program insert
```
./2opt <../Data/dataset.tsp> <number of climbers>
```
in terminal with a dataset and a given number og climbers.  

# Abstract
This project aims to efficiently solve the Travelling Salesman Problem (TSP) on the GPU. A short survey of three popular heuristic algorithms (2-Opt local search, Ant Colony Optimization, and Simulated annealing) is made. 2-Opt was chosen to be parallelized and implemented in CUDA. The 2 opt algorithm uses "climbers" that start with an initial solution, and improves it until a local optimum has been found. The GPU implementation exploits both degrees of parallelism by assigning one climber to a CUDA block (of threads). The code is publicly available at https://github.com/henni3/BA-project. The solution is compared with another 2-Opt implementation on the GPU by O’Neil et al. (https://userweb.cs.txstate.edu/~mb92/papers/pdpta11b.pdf) where one climber is assigned to one thread, i.e., this implementation uses only the parallelism of (different) climbers. Our project also evaluates how the accuracy of the solution depends on the number of different climbers. 
Compared to the implementation by O’Neil et al., this project achieved up to 27 times speed up on a 100-city data set with 1000 climbers but falls off for greater numbers of climbers. The accuracy of the solutions was also within < 5% of the optimal solution for 1000 climbers, indicating that a greater number of climbers for some data sets might not be necessary. Our implementation does not reach the big throughput of the GPU of 1555 GB/s for smaller data sets with < 100 cites, but for bigger data sets, the implementation exceeds the big throughput, since the majority of memory accesses happen to faster L2 cache. This shows that the implementation fully utilizes the GPU.

