#include <iostream>
#include <cmath>
#include <vector>
#include <complex>

void sieve();
std::vector<std::complex<double>> mixrad_fft(std::vector<double> &input);
std::vector<std::complex<double>> rad2_fft(std::vector<double> &input);
std::vector<double> periodogram(std::vector<std::complex<double>> &spectral);
