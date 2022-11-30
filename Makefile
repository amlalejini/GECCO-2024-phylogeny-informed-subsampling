# Project-specific settings
EMP_DIR := third-party/Empirical/include
CEREAL_DIR := third-party/Empirical/third-party/cereal/include

PROJECT := diagnostics
MAIN_CPP ?= source/diagnostics.cpp

# Flags to use regardless of compiler
CFLAGS_all := -Wall -Wno-unused-function -std=c++17 -lstdc++fs -I$(EMP_DIR)/ -Iinclude/ -Ithird-party/

# Native compiler information
CXX_nat := g++-12
CFLAGS_nat := -O3 -DNDEBUG -msse4.2 $(CFLAGS_all)
CFLAGS_nat_debug := -g $(CFLAGS_all)

default: $(PROJECT)
native: $(PROJECT)
all: $(PROJECT) $(PROJECT).js

debug:	CFLAGS_nat := $(CFLAGS_nat_debug)
debug:	$(PROJECT)

$(PROJECT): ${MAIN_CPP} include/
	$(CXX_nat) $(CFLAGS_nat) ${MAIN_CPP} -o $(PROJECT)

clean:
	rm -f $(PROJECT) web/$(PROJECT).js web/*.js.map web/*.js.map *~ source/*.o

# Debugging information
print-%: ; @echo '$(subst ','\'',$*=$($*))'
