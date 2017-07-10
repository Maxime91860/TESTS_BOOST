# BOOST_PATH=/ccc/products/boost-1.59.0/intel
# BOOST_PATH=/ccc/scratch/cont001/ocre/kermarquerm/boost_1_64_0/builded_libs
BOOST_PATH=/home/max/boost_1_64_0/install
BOOST_IFLAG=-I$(BOOST_PATH)/include
BOOST_LDFLAG=-L$(BOOST_PATH)/lib -lboost_serialization
BOOST_FLAGS=$(BOOST_IFLAG) $(BOOST_LDFLAG)

all : ex1.exe ex2.exe ex3.exe ex4.exe ex5.exe ex6.exe ex7.exe ex8.exe ex9.exe ex10.exe

%.exe : src/%.cpp
	c++ $^ $(BOOST_FLAGS) -g -std=c++11 -o bin/$@


clean : 
	rm -rf bin/* files_serialized/*
