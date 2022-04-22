COMPILER?=nvcc
OPT_FLAGS?=-O3 -arch=compute_35
LDFLAGS?=-lOpenCL

TEST_PBB=test-pbb
TEST_TSP=test-tsp
TEST_TSP2=test-tsp2

.PHONY: clean all run

default: compile run

compile: $(TEST_PBB) $(TEST_TSP) $(TEST_TSP2)

$(TEST_PBB): testPBB.cu hostSkel.cu.h pbbKernels.cu.h constants.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_PBB) testPBB.cu

$(TEST_TSP): tsp-kernels.cu.h tsp-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h transpose-host.cu.h transpose-kernels.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP) tsp-main.cu

$(TEST_TSP2):dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-kernels2.cu.h tsp-main2.cu
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP2) tsp-main2.cu

all: $(TEST_PBB) $(TEST_TSP) $(TEST_TSP2)

run-pbb: $(TEST_PBB)
	./$(TEST_PBB) 50003565 128

clean:
	rm -f $(TEST_PBB) $(TEST_TSP)
