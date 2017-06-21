BOOST_PATH=/ccc/products/boost-1.59.0/intel
# BOOST_PATH=/ccc/scratch/cont001/ocre/kermarquerm/boost_1_64_0/builded_libs
BOOST_IFLAG=-I$(BOOST_PATH)/include
BOOST_LDFLAG=-L$(BOOST_PATH)/lib -lboost_serialization
BOOST_FLAGS=$(BOOST_IFLAG) $(BOOST_LDFLAG)

all : ex1

ex1 : exemple1.cpp
	c++ $^ $(BOOST_FLAGS) -o $@.exe

clean : 
	rm -rf *.exe files_serialized/*