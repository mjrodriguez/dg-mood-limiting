BUILD_DIR=./build
INCLUDE_DIR=./include

CXX=g++
LFLAGS= -lblas -llapack
INCFLAGS= -I$(INCLUDE_DIR)

DEPS=$(BUILD_DIR)/.d
SOURCES=$(wildcard src/*.cpp)
OBJECTS=$(SOURCES:src/%.cpp=$(BUILD_DIR)/%.o)
EXEC=dg_mood_1d

.PHONY: all clean
all: $(EXEC)

# Build the executable
$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) $(LFLAGS) -o $@

$(BUILD_DIR)/%.o: Makefile | $(BUILD_DIR)
	$(CXX) -c $(INCFLAGS) -o $@ src/$*.cpp


$(DEPS): $(SOURCES) | $(BUILD_DIR)
	$(CXX) -c -std=c++11 $(INCFLAGS) -MM $(SOURCES) | awk '/^.*\.o:/{ print "$(BUILD_DIR)/" $$0; next } 1'  >$@

$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(EXEC) $(BUILD_DIR)