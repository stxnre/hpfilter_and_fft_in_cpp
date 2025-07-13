#include "fft.h"

const double pi = std::acos(-1);

void rad2_fft(std::vector<double> &input){
   /*
   Implementation of the Radix-2 Cooley-Tukey Fast Fourier Transform.
   The input size must be a positive power of 2.
   Can also perform the Inverse Fast Fourier Transform.
   */
   const int N = input.size();

   if(std::log2(N) % 1 != 0){
    cout << "Cannot Perform Radix-2 FFT with this input.";
    return;
   } 

   if(N==1){
    return input;
   }

   std::vector<std::complex<double>> out; //the output
   std::vector<double> odd;
   std::vector<double> even;
   int i,j;

   for(i=0;i<N/2;i+=2){
    odd[i] = input[i];
    even[i] = input[i+1];
   }

   std::vector<std::complex<double>> oddfft = rad2_fft(odd);
   std::vector<std::complex<double>> evenfft = rad2_fft(even);

   for(j = 0;j<N/2;j++){
        std::complex<double> t = std::polar(1.0, -2 * pi * j / N) * odd[j];
        out[j] = (even[j] + t);
        out[j + N/2] = (even[j] - t);
   }
   return out;
}


std::vector<double> periodogram(std::vector<std::complex<double>> &spectral){
    int i;
    const int N = spectral.size();
    vector<double> result(N);
    for(i=0;i<N;i++){
        result[i] = pow(abs(spectral[i]),2) / N;
    }
    return result;
}