COMPILER?=nvcc
OPT_FLAGS?=-O3
LDFLAGS?=-lOpenCL

TEST_PBB=test-pbb
TEST_TSP=test-tsp

.PHONY: clean all run

default: compile run

compile: $(TEST_PBB) $(TEST_TSP)

$(TEST_PBB): testPBB.cu hostSkel.cu.h pbbKernels.cu.h constants.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_PBB) testPBB.cu

$(TEST_TSP): tsp-kernels.cu.h tsp-main.cu hostSkel.cu.h pbbKernels.cu.h constants.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP) tsp-main.cu

all: $(TEST_PBB) $(TEST_TSP)

run-pbb: $(TEST_PBB)
	./$(TEST_PBB) 50003565 128

clean:
	rm -f $(TEST_PBB) $(TEST_TSP)
