# Time Series Seasonal-Trend Decomposition (in C++)

This repository contains code that performs Seasonal-Trend Decomposition on a time series, using the Hodrick-Prescott Filter and the Radix-2 Fast Fourier Transformation. I wrote these from scratch in C++ as an exercise. The implementations are contained within the folder `fft_hpfilter`, along with the executable (when compiled) that performs the time series decomposition.

The python scripts (to be written) will visualize the output of the C++ code.


## Part 1: Running C++ Code

Navigate to `fft_hpfilter`. It is recommended that you are using a virtual environment handler like mambas, as you can compile the C++ code just as I did. 

Once in the folder, type `make` to compile the file. `make clean` will remove the compiled file.

## Part 2: Python Visualizations

