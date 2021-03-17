include ../../modules/makefiles/definitions

INSTALL_DIR = $(TOOLS_DIR)/bin

# Additional fortran compilation flags
FFLAGS =  -O

# Additional linker flags
LFLAGS =  -lmodules

# Targets ...

TARGETS = pgmave pgm2ni ni2pgm autocorr raw2pgm pgm2mask circave braggextract

all: $(TARGETS)

braggextract: braggextract.o cluster_functions.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

circave: circave.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

pgm2mask: pgm2mask.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

pgmave: pgmave.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

pgm2ni: pgm2ni.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

ni2pgm: ni2pgm.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

raw2pgm: raw2pgm.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

autocorr: autocorr.o
	$(FC) $(FFLAGS) -o $@ $^ $(LFLAGS)

install: $(TARGETS)
	cp $(TARGETS) $(INSTALL_DIR)

clean:
	rm -f *.$(OBJSUFFIX) *.$(MODSUFFIX) $(TARGETS)
