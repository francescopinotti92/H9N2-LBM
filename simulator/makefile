# determine the platform
UNAME_S := $(shell uname -s)

# CXX and CXXFLAGS
ifeq ($(UNAME_S),Darwin)
  CXX := clang++
  CXXFLAGS:=-Ofast -Wall -shared -std=c++11 -undefined dynamic_lookup -fPIC
else
  CXX := g++
  CXXFLAGS:=-Ofast -Wall -shared -std=c++11 -undefined -fPIC
endif

# pybind11 include paths
INCPYBIND := $(shell python3 -m pybind11 --includes)

# all include paths
INC:=$(INCPYBIND) -I/src -I${CURDIR}/src

# use output of python3 -config --extension-suffix
#EXT:=.cpython-38-darwin.so
EXT:=$(shell python3-config --extension-suffix)

DEPS:=src/random.cpp src/simulator.cpp src/simulate_experiment.cpp src/bind_simulate_experiment.cpp

lib/simulate_experiment.$(EXT): $(DEPS)
	$(CXX) $(CXXFLAGS) $(INC) $(DEPS) -o lib/simulate_experiment$(EXT)

clean:
	@rm lib/*$(EXT)

.PHONY=clean

