CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -Wpedantic
INC = -Icpp/include
HDRS = $(wildcard cpp/include/radar/*.hpp)

all: build/radar_sim build/test_radar

build/radar_sim: cpp/apps/radar_sim.cpp $(HDRS) | build
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@

build/test_radar: cpp/tests/test_radar.cpp $(HDRS) | build
	$(CXX) $(CXXFLAGS) $(INC) $< -o $@

build:
	mkdir -p build

test: build/test_radar
	./build/test_radar

run: build/radar_sim
	./build/radar_sim

clean:
	rm -rf build

.PHONY: all test run clean
