MFEM_DIR = $(DIR_mfem)
MFEM_BUILD_DIR = $(MFEM_DIR)
CONFIG_MK = $(MFEM_BUILD_DIR)/config/config.mk

-include $(CONFIG_MK)

# Overall structure
BUILD_DIR=build
INCLUDE_DIR=./
DEPS=$(BUILD_DIR)/.d
SOURCES=$(wildcard impl/*.cpp)
OBJECTS=$(SOURCES:impl/%.cpp=$(BUILD_DIR)/%.o)
EXEC=dg_mood

CXXFLAGS=-g $(MFEM_CXXFLAGS) -DMFEM_DIR=\"$(MFEM_DIR)\"

# Use compiler configuration from MFEM
LFLGAS=${MFEM_LIBS} ${MFEM_EXT_LIBS} -lcublas -lblas -llapack -lgfortran
INCFLAGS=${MFEM_INCFLAGS} -I$(INCLUDE_DIR)

.PHONY: all clean style
all: $(EXEC)

# Build the executable
$(EXEC): $(OBJECTS)
	$(MFEM_CXX) $(OBJECTS) $(LFLGAS) -o $@

$(BUILD_DIR)/%.o: Makefile | $(BUILD_DIR)
	$(MFEM_CXX) -c $(CXXFLAGS) $(INCFLAGS) -o $@ impl/$*.cpp

# Use the compiler to determine dependencies on header files
# Some awk magic in the next target
# Prefix all lines matching the regexp /^.*\.o/ with build/ (to match e.g. "file.o:")
# and leave other lines alone
$(DEPS): $(SOURCES) | $(BUILD_DIR)
	$(MFEM_CXX) -c -std=c++11 $(INCFLAGS) -isystem $(MFEM_DIR) -MM $(SOURCES) \
		| awk '/^.*\.o:/{ print "$(BUILD_DIR)/" $$0; next } 1'  >$@

# Rebuild dependencies unless doing "make clean" or "make style"
ifneq ($(MAKECMDGOALS),clean)
ifneq ($(MAKECMDGOALS),style)
-include $(DEPS)
endif
endif

# Ensure build directory exists
$(BUILD_DIR):
	mkdir -p $@

clean:
	rm -rf $(EXEC) $(BUILD_DIR)

FORMAT_FILES = $(wildcard impl/*.?pp) $(wildcard ./*.?pp)
ASTYLE = astyle --options=$(MFEM_DIR)/config/mfem.astylerc
style:
	@if ! $(ASTYLE) $(FORMAT_FILES) | grep Formatted; then\
	   echo "No source files were changed.";\
	fi