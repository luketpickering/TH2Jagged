TH2JaggedTest: TH2JaggedTest.cxx TH2Jagged.h
	g++ $< -o $@ -g -O0 `root-config --cflags --glibs`
