CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall

SIMDIR   := simulations
V1_SRC   := $(SIMDIR)/FEA_sim_v1.cpp
V2_SRC   := $(SIMDIR)/FEA_sim_v2.cpp
V1_BIN   := FEA_sim_v1
V2_BIN   := FEA_sim_v2

.PHONY: all clean run run-v1 run-v2 v1 v2

all: $(V2_BIN)

v1: $(V1_BIN)
v2: $(V2_BIN)

$(V1_BIN): $(V1_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

$(V2_BIN): $(V2_SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

run: run-v2

run-v1: $(V1_BIN)
	./$(V1_BIN)

run-v2: $(V2_BIN)
	./$(V2_BIN)

clean:
	rm -f $(V1_BIN) $(V2_BIN)
