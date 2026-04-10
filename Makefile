CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall

TARGET   := fea_sim
SRC      := FEA_Chip_sim.cpp

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

clean:
	rm -f $(TARGET)
