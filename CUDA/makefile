COMPILER?=nvcc
OPT_FLAGS?=-O3 --use_fast_math
LDFLAGS?=-lOpenCL

TEST_TSP=test-tsp
REP_TSP=rep-tsp
.PHONY: clean all run

default: compile run

compile: $(TEST_TSP) $(REP_TSP) $(TEST_FUNC)

$(TEST_TSP): tsp-kernels.cu.h tsp-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP) tsp-main.cu

$(REP_TSP):dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-kernels.cu.h tsp-main-testing.cu tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(REP_TSP) tsp-main-testing.cu

$(TEST_FUNC):dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-kernels.cu.h main-test-func.cu tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_FUNC) main-test-func.cu

all: $(TEST_TSP) $(REP_TSP) $(TEST_FUNC)


clean:
	rm -f $(TEST_TSP) $(REP_TSP) $(TEST_FUNC)
