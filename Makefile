CXX = g++

# --- CONFIGURATION ARMADILLO ---
# On demande à pkg-config où se trouvent les headers (-I...)
ARMA_INC = $(shell pkg-config --cflags armadillo)
# On demande à pkg-config où se trouvent les librairies (-L... -larmadillo)
ARMA_LIBS = $(shell pkg-config --libs armadillo)

# Si pkg-config ne renvoie rien (cas rare), on garde une valeur par défaut
ifeq ($(ARMA_LIBS),)
    ARMA_LIBS = -larmadillo
endif

# --- FLAGS DE COMPILATION ---
# On ajoute $(ARMA_INC) ici pour que la compilation des .o trouve <armadillo>
CXXFLAGS = -Wall -O3 -std=c++14 $(CXXFLAGS_diff) $(ARMA_INC)
LDFLAGS = $(LDFLAGS_mac)

# On utilise les libs trouvées par pkg-config pour l'édition de lien
LIBS = $(ARMA_LIBS)

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

# Édition de liens (Création de l'exécutable)
$(TARGET): $(OBJS) | $(BIN_DIR)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(TARGET) $(LIBS)

# Compilation des objets (.cpp -> .o)
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
	rm -rf output/
	rm -f *.png
	# Si le dossier test existe et a un Makefile :
	# $(MAKE) -C test clean

# Compilation et exécution des tests unitaires
.PHONY: test
test:
	$(MAKE) -C test
	./test/test_poly

run: $(TARGET)
	./$(TARGET)