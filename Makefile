CXX      ?= c++
CXXFLAGS ?= -std=c++17 -O2 -Wall

TARGET   := FEA_sim
SRC      := FEA_sim.cpp

.PHONY: all clean run

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $@ $<

run: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
