CXX=g++
CXXTEST=cxxtestgen
CXXFLAGS= -Wall $(CXXFLAGS_diff)
LDFLAGS = $(LDFLAGS_mac)*
LIBS = -larmadillo

POLY_OBS=poly.o
TEST_CPP=t.cpp 
TEST_OBJ=t.o 

BIN_DIR=../bin 
SRC_DIR=../src 
TEST_SRC=test_poly.h
POLY_SRC=$(SRC_DIR)/poly.cpp
TARGET = test_poly

all:$(TARGET)

$(TEST_CPP):$(TEST_SRC)
	$(CXXTEST) --error-printer -o $(TEST_CPP) $(TEST_SRC)

$(POLY_OBJ):$(POLY_SRC)
	$(CXX) $(CXXFLAGS) -c $(POLY_SRC) -o $(POLY_OBJ)

$(TEST_OBJ):$(TEST_CPP)
	$(CXX) $(CXXFLAGS) -c $(TEST_CPP) -o $(TEST_OBJ)     

$(TARGET):$(TEST_OBJ) $(POLY_OBJ)
	$(CXX) $(LIBS) $(TEST_OBJ) $(POLY_OBJ) -o $(TARGET) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TEST_CPP)

