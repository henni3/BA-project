COMPILER?=nvcc
OPT_FLAGS?=-O3 --use_fast_math
LDFLAGS?=-lOpenCL

TSP=tsp
TEST=test
TEST_GB=test_gb
PRE_PRO=pre_pro

.PHONY: clean all run

default: compile run

compile: $(TSP) $(TEST) $(TEST_GB) $(PRE_PRO)

$(TSP): tsp-kernels.cu.h tsp-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TSP) tsp-main.cu

$(TEST): tsp-kernels.cu.h tsp-testing-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST) tsp-testing-main.cu
	
$(TEST_GB): tsp-kernels.cu.h tsp-testing-main.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_GB) tsp-gb-test.cu

$(PRE_PRO): tsp-kernels.cu.h tsp-run-program.cu dataCollector.cu.h hostSkel.cu.h pbbKernels.cu.h constants.cu.h tsp-main-helper.cu.h tsp-ker-helper.cu.h
	$(COMPILER) $(OPT_FLAGS) -o $(PRE_PRO) tsp-run-program.cu

all: $(TSP) $(TEST) $(TEST_GB) $(PRE_PRO)


clean:
	rm -f $(TSP) $(TEST) $(TEST_GB) $(PRE_PRO)
