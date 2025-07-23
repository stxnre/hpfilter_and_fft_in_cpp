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
    int lambda = 1000;
    sieve();

    // Read data
    std::ifstream inputFile("aus_production.csv");
    if(!inputFile.is_open()) throw std::runtime_error("Error: Could not open file.");

    std::string header;
    std::getline(inputFile, header);  // header not needed

    std::vector<std::string> dates;
    std::vector<std::vector<double>> data;

    std::string line;

    while (std::getline(inputFile, line)) {
        std::stringstream ss(line);
        std::string token;

        // Read date column
        std::getline(ss, token, ',');
        dates.push_back(token);

        std::vector<double> row;
        for (int i = 0; i < 6; ++i) {
            std::getline(ss, token, ',');
            
            if (token=="NA"||token.empty()) {
                row.push_back(std::nan(""));
            } else {
                row.push_back(std::stod(token));
            }
        }
        data.push_back(row);
    }

    inputFile.close();

    // Electricity Production
    std::vector<double> beer;
    for (const auto& row : data) {
    if (row.size() > 4) {
        beer.push_back(row[0]);
    }
}

    // Smoothing
    std::pair<std::vector<double>,std::vector<double>> smoothed = hpfilter_lapacke(beer,lambda);
    std::vector<double>* trend = &std::get<0>(smoothed);
    std::vector<double>* seasonal = &std::get<1>(smoothed);

    // Save seasonal-trend decomposition
    std::ofstream outputFile("season_trend_decomp.csv");
    if(!outputFile.is_open()) throw std::runtime_error("Error: Could not open output file.");

    outputFile << "Quarter,Original,Trend,Seasonal\n";
    for(int i = 0;i<dates.size();++i){
        outputFile << dates[i] << ','
                    << beer[i] << ','
                    << (*trend)[i] << ','
                    << (*seasonal)[i] << '\n';
    }
    outputFile.close();

    // Spectral Analysis
    std::vector<std::complex<double>> freq_rep = mixrad_fft(time_rep);
    std::vector<double> psd = periodogram(freq_rep); 

    // Save PSD result
    std::ofstream periodFile("periodogram.csv");
    if(!periodFile.is_open()) throw std::runtime_error("Error: Could not open output file.");

    for(int l=0;l<psd.size();++l){
        periodFile << psd[l] << '\n';
    }

    periodFile.close();

    return 0;
}