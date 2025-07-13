#include <cmath>
#include <vector>
#include <utility>
#include <string>
#include <iostream>


//function declarations
void pentadiag_cholesky_decomp(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c );
std::vector<double> pentadiag_cholesky_solver(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c, std::vector<double> &z );
std::pair<std::vector<double>,std::vector<double>> hpfilter(const std::vector<double> &series,const int lambda);