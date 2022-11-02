CXX = g++
CXXFLAGS = -O3 -Wall -Wextra -pedantic -std=c++14 -fopenmp

all: diffusion

diffusion: diffusion.cpp
	$(CXX) $< $(CXXFLAGS) -o $@

diffusion_naive: diffusion_naive.cpp
	$(CXX) $< $(CXXFLAGS) -o $@

run: diffusion
	./diffusion 100 1

jobs: diffusion
	./run.sh

stats: diffusion
	./stats.sh


clean:
	rm -rf diffusion lsf*

.PHONY: all clean
.SUFFIXES:
