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
            row.push_back(std::stod(token));
        }
        data.push_back(row);
    }

    inputFile.close();

    // Electricity Production
    std::vector<double> x;
    for (const auto& row : data) {
    if (row.size() > 4) {
        x.push_back(row[4]);
    }
}

    // Smoothing
    std::pair<std::vector<double>,std::vector<double>> smoothed = hpfilter_lapacke(x,lambda);
    std::vector<double>* trend = &std::get<0>(smoothed);
    std::vector<double>* seasonal = &std::get<1>(smoothed);

    // Save seasonal-trend decomposition
    std::ofstream outputFile("season_trend_decomp.csv");
    if(!outputFile.is_open()) throw std::runtime_error("Error: Could not open output file.");

    outputFile << "Quarter,Original,Trend,Seasonal\n";
    for(int i = 0;i<dates.size();i++){
        outputFile << dates[i] << ','
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