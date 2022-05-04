COMPILER?=nvcc
OPT_FLAGS?=-O3 -arch=compute_35 --use_fast_math
LDFLAGS?=-lOpenCL


TEST_PBB=test-pbb
TEST_TSP=test-tsp
TEST_TSP2=test-tsp2
TEST_TSP100=test-tsp100

.PHONY: clean all run

default: compile run

compile: $(TEST_PBB) $(TEST_TSP) $(TEST_TSP2) $(TEST_TSP100)

$(TEST_PBB): testPBB.cu hostSkel.cu.h pbbKernels.cu.h constants.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_PBB) testPBB.cu

$(TEST_TSP): tsp-kernels.cu.h tsp-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP) tsp-main.cu

$(TEST_TSP2):dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-kernels2.cu.h tsp-main2.cu
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP2) tsp-main2.cu

$(TEST_TSP100):dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-kernels-100.cu.h tsp-main-100.cu
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP100) tsp-main-100.cu

all: $(TEST_PBB) $(TEST_TSP) $(TEST_TSP2) $(TEST_TSP100)

run-pbb: $(TEST_PBB)
	./$(TEST_PBB) 50003565 128

clean:
	rm -f $(TEST_PBB) $(TEST_TSP)
