# ==========================================
# IPS-PROD Project Makefile
# ==========================================

# 1. Directories
SRC_DIR  = src
INC_DIR  = include
OBJ_DIR  = obj
TEST_DIR = test

# 2. Compiler Settings
CXX      = g++
# -std=c++11 : Mandatory language standard
# -Wall      : Turn on all standard warnings
# -Wextra    : Turn on extra warnings
# -O3        : Maximum optimization
# -I$(INC_DIR): Look for headers in 'include/'
CXXFLAGS = -std=c++11 -Wall -Wextra -O3 -I$(INC_DIR)

# 3. Libraries
LDFLAGS  = -larmadillo

# 4. Source Files
SRCS     = $(SRC_DIR)/main.cpp \
           $(SRC_DIR)/poly.cpp \
           $(SRC_DIR)/basis.cpp \
           $(SRC_DIR)/solver.cpp \
           $(SRC_DIR)/writer.cpp

# 5. Object Files
# Maps src/filename.cpp -> obj/filename.o
OBJS     = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.o, $(SRCS))

# 6. Targets
TARGET      = main
TEST_RUNNER = $(TEST_DIR)/runner.cpp
TEST_TARGET = $(TEST_DIR)/runner

# ==========================================
# Rules
# ==========================================

all: $(TARGET)

# Create the object directory if it doesn't exist
$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

# --- MAIN APPLICATION ---

# Link Rule
$(TARGET): $(OBJS)
	@echo "Linking $(TARGET)..."
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LDFLAGS)

# Compilation Rule
# The "| $(OBJ_DIR)" part is crucial: it ensures the directory exists before compiling
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) -c $< -o $@

# --- UNIT TESTS ---

# 1. Generate runner
$(TEST_RUNNER): $(TEST_DIR)/test_mandatory.hpp
	@echo "Generating test runner..."
	cxxtestgen --error-printer -o $@ $<

# 2. Compile and Run Tests
# Depends on $(OBJS) so it will trigger the compilation rule (and folder creation) above
tests: $(TEST_RUNNER) $(OBJS)
	@echo "Building tests..."
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_RUNNER) \
	$(filter-out $(OBJ_DIR)/main.o, $(OBJS)) $(LDFLAGS)
	@echo "Running tests..."
	@./$(TEST_TARGET)

# --- CLEANUP ---

clean:
	@echo "Cleaning up..."
	rm -rf $(OBJ_DIR) $(TARGET) $(TEST_TARGET) $(TEST_RUNNER)

.PHONY: all clean tests
