#include "hpfilter.h"

void pentadiag_cholesky_decomp(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c ){
    /*
    Derives the Cholesky Decomposition for a pentadiagonal symmetric matrix, represented by a,b,and c.
    Here, a is the diagonal of size n, b is the first upper diagonal of size n-1, and c is
    the second upper diagonal of size n-2.
    
    The vectors will be modified in place.
    */
    int N = a.size();

    // base case
    a[0] = std::sqrt(a[0]);
    b[0] /= a[0];
    a[1] = std::sqrt(a[1] - std::pow(b[0],2));

    for(int i=2;i<N;++i){
        c[i-2] /= a[i-2];
        b[i-1] = (b[i-1] - (b[i-2]* c[i-2])) / a[i-1];
        a[i] = std::sqrt(a[i] - std::pow(c[i-2],2) - std::pow(b[i-1],2));
    }
    return;
}

std::vector<double> pentadiag_cholesky_solver(std::vector<double> &a, std::vector<double> &b, std::vector<double> &c,const std::vector<double> &z ){
    /*
    Solves linear system Ax = z, where A is symmetric and pentadiagonal (and has been Cholesky-decomposed and stored in
    a,b, and c). 

    If run successfully, will output a new vector. 
    */
    int N = a.size();

    std::vector<double> x(N,0);

    //forward substitution
    x[0] = z[0] / a[0];
    x[1] = (z[1] - (b[0] *x[0])) / a[1];
    for(int i=2;i<N;++i){
        x[i] = (z[i] - (c[i-2]*x[i-2]) - (b[i-1]*x[i-1])) / a[i];
    }

    //backward substitution
    x[N-1] /= a[N-1];
    x[N-2] = (x[N-2] - (b[N-2]* x[N-1]))/ a[N-2];
    for(int i=N-3;i>=0;--i){
        x[i] = (x[i] - (c[i]*x[i+2]) - (b[i]*x[i+1])) / a[i];
    }

    return x;
}


std::pair<std::vector<double>,std::vector<double>> hpfilter(const std::vector<double> &series,const int lambda = 1600){
    /*
    The Hodrick-Prescott Filter to extract the trend from a time series.
    This function will provide both the trend and seasonal component.
    */
    
    int N = series.size();

    //second differences I + lambda(DtD)
    std::vector<double> d0(N,1 + 6 * lambda);
    std::vector<double> d1(N-1,-4 * lambda);
    std::vector<double> d2(N-2,lambda);

    d0[0] = 1 + lambda;
    d0[1] = 1 + 5 * lambda;
    d0[N-1] = 1 + 5 * lambda;
    d1[0] = -2 * lambda;

    pentadiag_cholesky_decomp(d0,d1,d2);

    std::vector<double> trend = pentadiag_cholesky_solver(d0,d1,d2,series);
    std::vector<double> seasonal(N);

    for(int i=0;i<N;++i)
    {
        seasonal[i] = series[i] - trend[i];
    }

    return make_pair(trend, seasonal);
}

std::pair<std::vector<double>,std::vector<double>> hpfilter_lapacke(const std::vector<double> &series, const int lambda = 1600){
    int kd = 2; // pentadiagonal system
    int Nrhs = 1; //series is one-dimensional
    int N = series.size();

    //second differences I + lambda(DtD)
    std::vector<double> d0(N,1 + 6 * lambda);
    std::vector<double> d1(N-1,-4 * lambda);
    std::vector<double> d2(N-2,lambda);

    d0[0] = 1 + lambda;
    d0[1] = 1 + 5 * lambda;
    d0[N-2] = 1 + 5 * lambda;
    d0[N-1] = 1 + lambda;
    d1[0] = -2 * lambda;
    d1[N-2] = -2 * lambda;


    //Banded matrix
    std::vector<double> banded((kd + 1) * N,0);
    for (int i = 0; i < N;++i) {
        banded[(kd) + i*(kd+1)] = d0[i];  // main diagonal on row kd
        if (i < N-1)
            banded[(kd-1) + (i+1)*(kd+1)] = d1[i]; // first superdiag on row kd-1, shifted right by 1
        if (i < N-2)
            banded[(kd-2) + (i+2)*(kd+1)] = d2[i]; // second superdiag on row kd-2, shifted right by 2
    }

    std::vector<double> series_copy = series;

    //Solve pentadiagonal system with LAPACKE
    int info = LAPACKE_dpbsv(LAPACK_COL_MAJOR, 'U', N, kd, Nrhs, banded.data(), kd + 1, series_copy.data(), N);

    if(info!=0){
        std::cout << "DPBSV failed with error code " << info << std::endl;
    }

    std::vector<double> seasonal(N);

    for(int i=0;i<N;++i)
    {
        seasonal[i] = series[i] - series_copy[i];
    }

    return make_pair(series_copy,seasonal);
}