COMPILER?=nvcc
OPT_FLAGS?=-O3 --use_fast_math
LDFLAGS?=-lOpenCL

TSP=tsp
TEST=test

.PHONY: clean all run

default: compile run

compile: $(TSP) $(TEST)

$(TSP): tsp-kernels.cu.h tsp-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TSP) tsp-main.cu

$(TEST): tsp-kernels.cu.h tsp-testing.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST) tsp-testing.cu

all: $(TSP) $(TEST)


clean:
	rm -f $(TSP) $(TEST)
