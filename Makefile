CXX = g++
CXXFLAGS = -Wall -O3 -std=c++14 $(CXXFLAGS_diff)
LDFLAGS = $(LDFLAGS_mac)
LIBS = -larmadillo
SRC_DIR = src
OBJ_DIR = obj
BIN_DIR = bin
INC_DIR = include

TARGET = $(BIN_DIR)/main

# Recherche tous les fichiers .cpp dans le dossier SRC_DIR
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
# Transforme tous les chemins de $(SRC_DIR) en chemin .o dans le dossier $(OBJ_DIR)
OBJS = $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRCS))


all: $(TARGET)


$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp | $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -I$(INC_DIR) -c $< -o $@

# Créer un dossier bin s'il n'existe pas
$(BIN_DIR):
	mkdir -p $(BIN_DIR)

# Créer un dossier obj s'il n'existe pas
$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

# Nettoyage
.PHONY: clean
clean:
	rm -f $(OBJ_DIR)/*.o
	rm -f $(TARGET)
	$(MAKE) -C test clean

# Compilation et éxécution des tests unitaires
.PHONY: tests
tests:
	$(MAKE) -C tests
	./tests/test_poly


run:$(TARGET)
	./$(TARGET)



