
BOOST_PATH=/ccc/products/boost-1.59.0/intel
BOOST_IFLAG=-I$(BOOST_PATH)/include
BOOST_LDFLAG=-L$(BOOST_PATH)/lib -lboost_serialization
BOOST_FLAGS=$(BOOST_IFLAG) $(BOOST_LDFLAG)

all : note.exe gps.exe demo.exe  test_perso.exe # demo_xml.exe

%.exe : %.cpp
	c++ $^ $(BOOST_FLAGS) -o $@

clean :
	rm -rf *.exe *.ser *.txt