MScG15-MatlabMex
======

Mex-Files for MScG15 matlab_tools

## Functions:
* [legendreFunctions](#legendreFunctions)

#### legendreFunctions
Computation of the Legendre functions and arranging them in a lower 
triangular-matrix. Speeds up the whole SH synthesis by ~25%.
* legendreFunctions.m
	* m-function wrapper for parsing inputs
* legendreFunctionsmx.c
	* c-implementation for computing the Legendre functions

## Installation
#### Requirements:
* [Matlab Compatible Compiler](https://de.mathworks.com/support/compilers.html)
* [GNU Scientific Library](https://www.gnu.org/software/gsl/)

#### Building MEX file:

Build MEX file by calling Matlabs *mex*-command and linking needed librarys. 

#### Replace original Matlab function
Copy the binary MEX file and m-function wrapper file to the
matlab_tools-directory. MEX files act as drop-in replacement for the original
functions.
