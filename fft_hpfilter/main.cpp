#include "hpfilter.h"
#include "fft.h"
#include <string>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <sstream>

/*

g++ main.o hpfilter.o fft.o -o seasonal_trend

*/

int main(){
    int lambda = 2000;

    // Read data
    std::ifstream inputFile("AEP_hourly.csv");
    if(!inputFile.is_open()) throw std::runtime_error("Error: Could not open file.");

    std::string header;
    std::getline(inputFile, header);  // header not needed

    std::vector<std::string> datetimes;
    std::vector<double> x;

    std::string line;
    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string datetime, x_str;

        if (std::getline(ss, datetime, ',') && std::getline(ss, x_str)) {
            datetimes.push_back(datetime);
            x.push_back(std::stod(x_str));
        }
    }
    inputFile.close();

    // Smoothing
    std::pair<std::vector<double>,std::vector<double>> smoothed = hpfilter_lapacke(x,lambda);
    std::vector<double>* trend = &std::get<0>(smoothed);
    std::vector<double>* seasonal = &std::get<1>(smoothed);

    // Save seasonal-trend decomposition
    std::ofstream outputFile("season_trend_decomp.csv");
    if(!outputFile.is_open()) throw std::runtime_error("Error: Could not open output file.");

    outputFile << "datetime,original,trend,seasonal\n";
    for(int i = 0;i<datetimes.size();i++){
        outputFile << datetimes[i] << ','
                    << x[i] << ','
                    << (*trend)[i] << ','
                    << (*seasonal)[i] << '\n';
    }
    outputFile.close();

    // Spectral Analysis
    int k = (int)std::log2(x.size());
    int N = pow(2,k);

    std::vector<double> time_rep(N);
    for(int j=0;j<N;j++){
        time_rep[j] = (*seasonal)[j];
    }

    std::vector<std::complex<double>> freq_rep = rad2_fft(time_rep);
    std::vector<double> psd = periodogram(freq_rep); 

    // Save PSD result
    std::ofstream periodFile("periodogram.csv");
    if(!periodFile.is_open()) throw std::runtime_error("Error: Could not open output file.");

    for(int l=0;l<psd.size();l++){
        periodFile << psd[l] << '\n';
    }

    periodFile.close();

    return 0;
}