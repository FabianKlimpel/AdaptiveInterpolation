# AdaptiveInterpolation
First rough implementation of an adaptive polynomial interpolation for MC generator tuning based on  
B. Zheng, J. Takamatsu, K. Ikeuchi et al., An Adaptive and Stable Method for
Fitting Implicit Polynomial Curves and Surfaces, IEEE vol. 32 no. 3 (2010)
561-568.

The input requires mutliple files containing observables in a YODA format and returns a single file with the interpolation for each data point. The I/O is written such that it fits the Professor (arXiv:0907.2973v1) framework in the version 2.1.4.
