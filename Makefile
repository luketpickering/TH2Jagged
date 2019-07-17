lib/libTH2Jagged.so: src/TH2Jagged.h src/TH2Jagged.cxx TH2JaggedDict.cxx lib TH2JaggedDict_rdict.pcm
	g++ TH2JaggedDict.cxx -c -o TH2JaggedDict.o -Isrc -fPIC `root-config --cflags`
	g++ src/TH2Jagged.cxx -c -o TH2Jagged.o -Isrc -fPIC `root-config --cflags`
	g++ -shared -o $@ TH2JaggedDict.o TH2Jagged.o `root-config --libs`

.PHONY: clean

lib:
	-mkdir -p lib

bin:
	-mkdir -p bin

TH2JaggedDict.cxx TH2JaggedDict_rdict.pcm: src/TH2Jagged.h src/TH2JaggedLinkdef.h
	cd src; rootcling -f ../TH2JaggedDict.cxx -c TH2Jagged.h TH2JaggedLinkdef.h

bin/TH2JaggedBinningTest: test/TH2JaggedBinningTest.cxx src/TH2Jagged.h lib/libTH2Jagged.so bin
	g++ $< -o $@ -g -O0 -Isrc -Llib -lTH2Jagged `root-config --cflags --glibs`

bin/TH2JaggedErrorTest: test/TH2JaggedErrorTest.cxx src/TH2Jagged.h lib/libTH2Jagged.so bin
	g++ $< -o $@ -g -O0 -Isrc -Llib -lTH2Jagged `root-config --cflags --glibs`


bin/TH2JaggedTest: test/TH2JaggedTest.cxx src/TH2Jagged.h lib/libTH2Jagged.so bin
	g++ $< -o $@ -g -O0 -Isrc -Llib -lTH2Jagged `root-config --cflags --glibs`

bin/TH2JaggedWriteTest: test/TH2JaggedWriteTest.cxx src/TH2Jagged.h lib/libTH2Jagged.so bin
	g++ $< -o $@ -g -O0 -Isrc -Llib -lTH2Jagged `root-config --cflags --glibs`

bin/TH2JaggedReadTest: test/TH2JaggedReadTest.cxx src/TH2Jagged.h lib/libTH2Jagged.so bin
	g++ $< -o $@ -g -O0 -Isrc -Llib -lTH2Jagged `root-config --cflags --glibs`

all: lib/libTH2Jagged.so bin/TH2JaggedTest bin/TH2JaggedWriteTest bin/TH2JaggedReadTest bin/TH2JaggedBinningTest bin/TH2JaggedErrorTest

clean:
	-rm -rf TH2JaggedDict.cxx bin lib *.pcm *.o *.so
