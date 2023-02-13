COMPILER?=nvcc
OPT_FLAGS?=-O3 --use_fast_math
LDFLAGS?=-lOpenCL

TEST_TSP=test-tsp
TEST_TSP100=test-tsp100
.PHONY: clean all run

default: compile run

compile: $(TEST_TSP) $(TEST_TSP100)

$(TEST_TSP): tsp-kernels.cu.h tsp-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP) tsp-main.cu

$(TEST_TSP100):dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-kernels.cu.h tsp-main-100.cu tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP100) tsp-main-100.cu

all: $(TEST_TSP) $(TEST_TSP100)


clean:
	rm -f $(TEST_TSP) $(TEST_TSP100)
