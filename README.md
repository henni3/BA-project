# Bachelor Project DIKU 2023 - Traveling Salesmens Problem
Supervisor: Associate Professor, Cosmin Eugen Oance  
Students: Henriette Naledi Winther Hansen and Jóhann Utne

# Structure of Github Repository
- The _report.pdf_ contains the decribtions of the project and related works, how to convert the sequential code to parallel, teori behind a GPU, a testing section and a conclussion based on the tests.  
- The _CUDA directory_ contains the CUDA implementation of the 2-Opt algorithm, together with publicly available code from O’Neil et.al (https://userweb.cs.txstate.edu/~mb92/papers/pdpta11b.pdf). It also contains the performance data from the performance tests on our version.  
- The _Data directory_, contains all the data sets acquired from the TSBLIB95 library (http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/tsp/).  
- The _Futhark directory_ contains the Futhark implementation of the 2-Opt algorithm.  
- The _Seqventielt directory_, contains the sequential implementations of the 2-Opt, ACO, and SA algorithms.  

# How to run the project
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

To compile the _sequential implementation_ go to the insert
```
make all
```
