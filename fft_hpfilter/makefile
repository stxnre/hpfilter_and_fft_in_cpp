FFT_LSTDFLG = -llapacke -llapack -lblas -lm
FFT_INCLUDE = -std=c++17 -I$(CONDA_PREFIX)/include
FFT_LIB = -L$(CONDA_PREFIX)/lib
FFT_OBJS = seasonal_decompose

all: ${FFT_OBJS}
	rm -f *.o

fft.o: fft.cpp fft.h
	clang++ -g -c fft.cpp -o fft.o ${FFT_INCLUDE}

hpfilter.o: hpfilter.cpp hpfilter.h
	clang++ -g -c hpfilter.cpp -o hpfilter.o ${FFT_INCLUDE}

main.o: main.cpp
	clang++ -g -c main.cpp -o main.o ${FFT_INCLUDE}

seasonal_decompose: main.o fft.o hpfilter.o
	clang++ main.o fft.o hpfilter.o -o seasonal_decompose ${FFT_LIB} ${FFT_LSTDFLG}

clean:
	rm -f *.o
	rm -f ${FFT_OBJS}