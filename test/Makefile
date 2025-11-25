CXX=g++
CXXTEST=cxxtestgen
CXXFLAGS= -Wall $(CXXFLAGS_diff)
LDFLAGS = $(LDFLAGS_mac)
LIBS = -larmadillo


POLY_OBJ=$(OBJ_DIR)/poly.o
TEST_CPP=test.cpp 
TEST_OBJ=test.o 

OBJ_DIR=../obj
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
	$(CXX)  $(TEST_OBJ) $(POLY_OBJ) -o $(TARGET) $(LIBS)

.PHONY: clean
clean:
	rm -f $(TEST_CPP) $(TEST_OBJ) $(POLY_OBJ)
	rm -f $(TARGET)



