all: clean larch

larch:
	g++ -O3 -std=c++17 -o larch_scan larch_scan.cpp  -lm

clean:
	/bin/rm -f larch_scan
