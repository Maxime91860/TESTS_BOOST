
BOOST_PATH=/home/max/boost_1_64_0/install
BOOST_IFLAG=-I$(BOOST_PATH)/include
BOOST_LDFLAG=-L$(BOOST_PATH)/lib -lboost_serialization
BOOST_FLAGS=$(BOOST_IFLAG) $(BOOST_LDFLAG)

all : note.exe gps.exe demo.exe  test_perso.exe # demo_xml.exe

%.exe : src/%.cpp
	c++ $^ $(BOOST_FLAGS) -o bin/$@

clean :
	rm -rf *.exe *.ser *.txt bin/*