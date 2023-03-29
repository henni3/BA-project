COMPILER?=nvcc
OPT_FLAGS?=-O3
LDFLAGS?=-lOpenCL

TEST_TSP=test-tsp11

.PHONY: clean all run

default: compile run

compile: $(TEST_TSP)


$(TEST_TSP): TSP-GPU11.cu
	$(COMPILER) $(OPT_FLAGS) -o $(TEST_TSP) TSP-GPU11.cu

all: $(TEST_TSP)

clean:
	rm -f $(TEST_TSP)
