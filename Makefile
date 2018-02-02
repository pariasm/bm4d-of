CFLAGS   = -std=c99 -O3 -DNDEBUG -w -fPIC 
CXXFLAGS = -fopenmp -O3 -DNDEBUG -w -fPIC 
CPPFLAGS = -std=c++11 -Ilib/iio -Isrc/core -Isrc/utils
LDLIBS   = -fopenmp -lcblas -ljpeg -lpng -ltiff -lm

BINDIR   = build/bin/
BIN      = $(BINDIR)nldct
SCRIPT   = $(BINDIR)nldct-mp4.sh
OBJ      = lib/iio/iio.o src/core/matrix_funs.o src/core/nldct.o \
			  src/utils/lib_image.o src/utils/lib_videot.o \
			  src/utils/mt19937ar.o src/utils/utilities.o src/main_vnlb.o

all      : $(BIN) $(SCRIPT)
$(BIN)   : $(OBJ) $(BINDIR) ; $(CXX) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)
$(SCRIPT): src/scripts/nldct-mp4.sh $(BINDIR); cp $< $@
$(BINDIR): ; mkdir -p $@
clean    : ; $(RM) $(BIN) $(OBJ) $(SCRIPT)

