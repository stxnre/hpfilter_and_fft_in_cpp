#include "fft.h"

const double pi = std::acos(-1);

#define MAXN = 100001;
vector<int> spfs(MAXN+1,1);


void sieve(){
    /*
    Executes Eratosthenes's Sieve to get the smallest prime factor for each number up to MAXN.
    */
    spf[0] = 0;
    for(int i = 2; i<=MAXN;++i){
        if(spf[i] == 1){
            for(int j = i;j<=MAXN;j+=i){
                if(spf[j]==1){
                    spf[j] = i;
                }
            }
        }
    }

}

std::vector<std::complex<double>> mixrad_fft(std::vector<double> &input){
    /*
    Implements Mixed-Radix Cooley-Tukey Fast Fourier Transform. 
    */

    const int N = input.size();
    std::vector<std::complex<double>> out(N);

    if(N==1){
        out[0] = std::complex<double>(input[0],0);
        return out;
    }

    const int radix = spf[N]; //smallest prime factor of N
    const int M = N / radix;
    std::vector<std::vector<std::complex<double>> X(radix,M);

    for(int r=0;r<radix;++r){

        std::vector<double> input_r(M);
        for(int k=0;k<M;++k){
            input_r[k] = input[r+(k*radix)];
        }

        std::vector<std::complex<double>> r_fft = mixrad_fft(input_r);
        for(int m=0;m<M;++m){
            std::complex<double> t = std::polar(1.0, -2 * pi * m * r / N) ;
            X[r][m] = t * r_fft[r];
        }
    }

    for(int r = 0;r<radix;++r){
        std::vector<std::complex<double>> temp(M, 0.0);
        for (int r = 0; r < radix; ++r) {
        std::complex<double> twiddle = std::polar(1.0,-2 * pi * k * M * r / N);
            for (int m = 0; m < M; ++m) {
            temp[m] += twiddle * X[r][m];
            }
        }
    // Copy temp into the correct slice of Output
        for (int m = 0; m < M; ++m) {
            out[k * M + m] = temp[m];
        }
    }
    return out;
}


std::vector<std::complex<double>> rad2_fft(std::vector<double> &input){
   /*
   Implementation of the Radix-2 Cooley-Tukey Fast Fourier Transform.
   The input size must be a non-negative power of 2.
   Can also perform the Inverse Fast Fourier Transform.
   */
   const int N = input.size();
   std::vector<std::complex<double>> out(N); //the output
   bool powerOfTwo = !(N == 0) && !(N & (N - 1));

   if(!powerOfTwo){
    std::cout << "Cannot Perform Radix-2 FFT with this input.";
    return out;
   } 

   if(N==1){
    out[0] = std::complex<double>(input[0],0);
    return out;
   }

   std::vector<double> odd(N/2);
   std::vector<double> even(N/2);

   for(int i=0;i<N/2;i++){
    odd[i] = input[2*i+1];
    even[i] = input[2*i];
   }

   std::vector<std::complex<double>> oddfft = rad2_fft(odd);
   std::vector<std::complex<double>> evenfft = rad2_fft(even);

   for(int j = 0;j<N/2;j++){
        std::complex<double> t = std::polar(1.0, -2 * pi * j / N) * odd[j];
        out[j] = (even[j] + t);
        out[j + N/2] = (even[j] - t);
   }
   return out;
}


std::vector<double> periodogram(std::vector<std::complex<double>> &spectral){
    const int N = spectral.size();
    std::vector<double> result(N);
    for(int i=0;i<N;i++){
        result[i] = pow(abs(spectral[i]),2) / N;
    }
    return result;
}