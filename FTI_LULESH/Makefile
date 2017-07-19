
#default build suggestion of MPI + OPENMP with gcc on Livermore machines you might have to change the compiler name

SHELL = /bin/sh
.SUFFIXES: .cc .o

LULESH_EXEC = lulesh2.0

MPI_INC = /opt/local/include/openmpi
MPI_LIB = /opt/local/lib

SERCXX = c++ -DUSE_MPI=0
MPICXX = mpic++ -DUSE_MPI=1
CXX = $(MPICXX)

SOURCES2.0 = \
	lulesh.cc \
	lulesh-comm.cc \
	lulesh-viz.cc \
	lulesh-util.cc \
	lulesh-init.cc
OBJECTS2.0 = $(SOURCES2.0:.cc=.o)
#:


#If you do not have the libraries used at the build you can get them at:
#silo:  https://wci.llnl.gov/codes/silo/downloads.html
#visit: https://wci.llnl.gov/codes/visit/download.html
#hdf5:  https://support.hdfgroup.org/downloads/
#boost: http://www.boost.org/users/download/
#fti:   https://github.com/leobago/fti

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#

#Change the paths in *_PATH variables

#Serialization with Boost library
BOOST_PATH=$(HOME)/boost_1_64_0/install
BOOST_IFLAG=-I$(BOOST_PATH)/include
BOOST_LDFLAG=-L$(BOOST_PATH)/lib -lboost_serialization
BOOST_FLAGS=$(BOOST_IFLAG) $(BOOST_LDFLAG)

#Checkpoint - Restart with FTI libray
FTIPATH=$(HOME)/FTI/install
FTI_IFLAG=-I$(FTIPATH)/include 
FTI_LDFLAG=$(FTIPATH)/lib/libfti.a

#HDF5 and Silo library to the vizulization part
# SILO_PATH=$(HOME)/silo-4.10.2
# SILO_IFLAG=-I$(SILO_PATH)/include
# SILO_LDFLAG=-Wl,-rpath -L$(SILO_PATH)/lib -Wl,$(SILO_PATH)/lib -lsiloh5
# HDF5_LDFLAG=-L$(HOME)/hdf5-1.8.19/lib -lhdf5

#----------------------------------------------------------------------------------------------#
#----------------------------------------------------------------------------------------------#


#Default build suggestions with OpenMP with checkpoint/restart
CXXFLAGS = -g -O3 -fopenmp -I. -Wall $(BOOST_IFLAG) $(FTI_IFLAG)
LDFLAGS = -g -O3 -fopenmp $(BOOST_LDFLAG) $(FTI_LDFLAG)

#Flags for using checkpoint/restart and to get vizulization
# CXXFLAGS = -g -O3 -fopenmp -DVIZ_MESH $(SILO_IFLAG) $(BOOST_IFLAG) $(FTI_IFLAG) -Wall -Wno-pragmas 
# LDFLAGS = -g -O3 -fopenmp $(SILO_LDFLAG) $(HDF5_LDFLAG) $(BOOST_LDFLAG) $(FTI_LDFLAG)

.cc.o: lulesh.h
	@echo "Building $<"
	$(CXX) -c $(CXXFLAGS) -o $@  $<

all: $(LULESH_EXEC)

lulesh2.0: $(OBJECTS2.0)
	@echo "Linking"
	$(CXX) $(OBJECTS2.0) $(LDFLAGS) -lm -o $@

clean:
	/bin/rm -f *.o *~ $(OBJECTS) $(LULESH_EXEC)
	/bin/rm -rf *.dSYM lulesh#* 
	/bin/rm -rf buffer_size#*
	make -C ckpt_files

tar: clean
	cd .. ; tar cvf lulesh-2.0.tar LULESH-2.0 ; mv lulesh-2.0.tar LULESH-2.0


run : lulesh2.0
	mpirun -n 8 ./lulesh2.0 -s 38 -p -i 1000
