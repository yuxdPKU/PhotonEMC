##############################################
# please add new classes in alphabetical order

AUTOMAKE_OPTIONS = foreign

AM_CPPFLAGS = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include

AM_LDFLAGS = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  -L$(OFFLINE_MAIN)/lib64 \
  -lSubsysReco \
  -lcalo_io \
  -lglobalvertex_io \
  -lphparameter \
  -lHepMC \
  -lCLHEP \
  -lcalo_io \
  -lg4detectors \
  -lg4detectors_io \
  -lmbd_io \
  -lglobalvertex \
  -lg4dst \
  -lg4eval

pkginclude_HEADERS = \
  PhotonEMC.h

# List of shared libraries to produce
lib_LTLIBRARIES = \
  libphotonemc.la

libphotonemc_la_SOURCES = \
  PhotonEMC.cc

libphotonemc_la_LIBADD = \
  -lphool \
  -lcalo_io \
  -lfun4all \
  -lg4detectors \
  -lg4detectors_io \
  -lcdbobjects \
  -lSubsysReco

################################################
# linking tests

noinst_PROGRAMS = \
  testexternals_photonemc

BUILT_SOURCES  = testexternals.cc

testexternals_photonemc_SOURCES = testexternals.cc
testexternals_photonemc_LDADD = libphotonemc.la

testexternals.cc:
	echo "//*** this is a generated file. Do not commit, do not edit" > $@
	echo "int main()" >> $@
	echo "{" >> $@
	echo "  return 0;" >> $@
	echo "}" >> $@

##############################################
# please add new classes in alphabetical order

clean-local:
	rm -f $(BUILT_SOURCES)
