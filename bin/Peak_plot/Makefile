OBJ_DIR = ./
BIN_DIR = ./
TABIX_DIR = ./
SOURCES= main.cpp
OBJECTS= main.o
BUILT_OBJECTS= $(patsubst %,$(OBJ_DIR)/%,$(OBJECTS))
PROGRAM= peak_plot
BIN_OBJECT   = $(patsubst %,$(BIN_DIR)/%,$(PROGRAM))

TABIX = tabix-0.2.6
INCLUDES= -I $(patsubst  %,$(TABIX_DIR)/%,$(TABIX))
LIBRARY = -L $(patsubst  %,$(TABIX_DIR)/%,$(TABIX))




LDFLAGS = -ltabix -lm -lz -lpthread

$(BIN_OBJECT): $(BUILT_OBJECTS)
	$(CXX) -o $@ $< $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)
$(BUILT_OBJECTS): $(SOURCES)
	@echo "  * compiling" $(*F).cpp
	@$(CXX)      -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)
	@echo $(CXX) -c -o $@ $(*F).cpp $(LDFLAGS) $(CXXFLAGS) $(DFLAGS) $(INCLUDES) $(LIBRARY)

clean:
	@echo "Cleaning up."
	@rm -f $(OBJ_DIR)/main.o $(BIN_OBJECT)

.PHONY: clean
