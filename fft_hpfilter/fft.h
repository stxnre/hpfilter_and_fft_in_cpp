#include <iostream>
#include <cmath>
#include <vector>
#include <complex>

std::vector<std::complex<double>> rad2_fft(std::vector<double> &input);
std::vector<double> periodogram(std::vector<std::complex<double>> &spectral);
