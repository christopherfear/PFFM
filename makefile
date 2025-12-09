CXX = g++
CXXFLAGS = -O3 -fopenmp -fno-math-errno -DEIGEN_NO_DEBUG -DNDEBUG -march=native \
           -I ~/PFFM/eigen-3.4.0 \
           -I /opt/SuiteSparse/include/suitesparse \
           -I /opt/OpenBLAS/include \
		   -I /usr/local/cuda/include
LDFLAGS = -L/opt/SuiteSparse/lib -L/opt/OpenBLAS/lib -Wl,-rpath,/opt/SuiteSparse/lib
LIBS = -lcholmod -lamd -lcolamd -lcamd -lccolamd -lsuitesparseconfig -lopenblas
TARGETI = Dcb2D-ModeI
TARGETII = Dcb2D-ModeII
SRCI = Dcb2D-ModeI.cpp
SRCII = Dcb2D-ModeII.cpp

all: $(TARGETI) $(TARGETII)

modeI: $(TARGETI)

modeII: $(TARGETII)

$(TARGETI): $(SRCI)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)

$(TARGETII): $(SRCII)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LIBS)
	
clean:
	-rm Dcb2D-ModeI Dcb2D-ModeII