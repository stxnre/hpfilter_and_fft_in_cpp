all: build

build:
	clang++ -std=c++17 main.cpp fft.cpp hpfilter.cpp -o seasonal_decompose -llapacke -llapack -lblas

clean:
	rm -f *.o seasonal_decompose