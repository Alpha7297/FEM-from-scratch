CXX=g++
CXXFLAGS=-O3 -Wall -std=c++11 -fopenmp
DEMO_TARGETS=demo
DEMO_SRC=gen.cpp LU.cpp cal.cpp demo.cpp
NEMO_TARGETS=nemo
NEMO_SRC=gen.cpp LU.cpp cal.cpp nemo.cpp
TEST_TARGET=test
TEST_SRC=gen.cpp test.cpp
DATA_FILES=points.txt triangles.txt phi.txt
.PHONY: demo nemo test clean

all: demo nemo

demo:$(DEMO_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^
	./demo
	python plotphi.py
nemo:$(NEMO_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^
	./nemo
	python plotphi.py
test:$(TEST_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $^
	./test
	python plotedge.py
clean:
	rm -f $(DEMO_TARGETS) $(DATA_FILES) $(TEST_TARGET) $(NEMO_TARGETS)