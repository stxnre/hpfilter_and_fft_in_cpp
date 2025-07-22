# Time Series Seasonal-Trend Decomposition (in C++)

This repository contains code that performs Seasonal-Trend Decomposition on a time series, using the Hodrick-Prescott Filter and the Radix-2 Fast Fourier Transformation. I wrote these from scratch in C++ as an exercise. The implementations are contained within the folder `fft_hpfilter`, along with the executable (when compiled) that performs the time series decomposition.

The python scripts (to be written) will visualize the output of the C++ code.

## The Data

`aus_production.csv` is taken from the `tsibbledata` package in R. It comprises of quarterly estimates of select indicators of manufacturing production in Australia, ranging from Q1 of 1956 to Q2 of 2010. 

The Hodrick-Prescott Filter is a method to smooth a time series, in order to try and decompose the trend and seasonal components. It is the solution to the following equation:

$$
\argmin_{\theta_t} \sum_{t=1}^n (x_t - \theta_t)^2 + \lambda \sum_{t=2}^{n-1}\left((\theta_{t+1} - \theta_t)- (\theta_t - \theta_{t-1})\right)^2
$$


## Part 1: Running C++ Code

Navigate to `fft_hpfilter`. It is recommended that you are using a virtual environment handler like mambas, as you can compile the C++ code just as I did. 

Once in the folder, type `make` to compile the file. `make clean` will remove the compiled file.

## Part 2: Python Visualizations

