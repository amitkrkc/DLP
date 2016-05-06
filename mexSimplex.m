% Description
% MEXSIMPLEX is a MEX file that tests the simplex.h file. Specifically, it solves the program y=argmin_{ y \in \Delta } 0.5*c*y^2-v^\top y where the scalar c can be negative or positive.
% 
% Compilation
% Compile it under MATLAB command as
% mex mexSimplex.cpp
% 
% Usage
% y=mexSimplex(v,c)
% where
%
% Input arguments
% 
% v is a n-by-1 dimensional vector. It must be a double.
% c is a scalar. It must be double.
%
%
% Output(s) 
% 
% y is a n-by-1 dimensional vector. It is returned as double.
%
%
% Written by Amit Kumar K C, PhD at Wollman Lab (UCSD)
% Date: 5th May 2016
