CXX = g++
COMMON_FLAGS =
CXXFLAGS = -Wall -msse3 -fexceptions -g -std=c++11 $(COMMON_FLAGS)
LINKER_FLAGS = $(COMMON_FLAGS)

SRCDIR = src
OUTDIR = bin
PRECOMPILED_HEADER = precompiled.h

CXXFLAGS += -include $(SRCDIR)/$(PRECOMPILED_HEADER)

EXEC = $(OUTDIR)/prg
SOURCES = $(wildcard src/*.cpp)
OBJECTS := $(SOURCES:$(SRCDIR)/%.cpp=$(OUTDIR)/%.o)
PRECOMPILED = $(SRCDIR)/$(PRECOMPILED_HEADER).gch

all: $(PRECOMPILED) $(EXEC)

$(PRECOMPILED): $(SRCDIR)/$(PRECOMPILED_HEADER)
	$(CXX) $(CXXFLAGS) $(SRCDIR)/$(PRECOMPILED_HEADER)

$(EXEC): $(OBJECTS)
	$(CXX) $(LINKER_FLAGS) $(OBJECTS) -o $(EXEC)

$(OBJECTS): $(OUTDIR)/%.o : $(SRCDIR)/%.cpp
	mkdir -p $(OUTDIR)
	$(CXX) $(CXXFLAGS) -H -c $< -o $@

clean:
	rm -rf $(OUTDIR)

.PHONY: clean
