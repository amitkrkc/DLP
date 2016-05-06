% Description
% MEXDLP is a MEX file that tests the dlp.h file. It performs a series of sanity checks before finally creating a DLP object. 
% 
% 
% Usage
% [Y, fevo]=mexGraph(parameters, edgelist, nucleilist, scorelist, lutlist)
% where
%
% Input arguments
% 
% parameters is a 1-by-11 dimensional (row) vector. It must be a double.
% edgelist is a m-by-3 dimensional matrix, where m is the number of edges. Each row consits of (source, target, weight). It must be double.
% nucleilist in a k-by-1 column vector, where k is the number of nuclei in the image. It must be int32.
% 
% scorelist is a n-by-1 column vector where n is the number of nodes. It must be double.
% 
% lutlist is a p-by-3 matrix where p is the number of superpixels that belong to the red cells. Each row consists of (superpixel index, red cell label, pixel area)
% 
% Output arguments
%
% Y is a n-by-k label assignent matrix where n is the number of nodes, and k is the label vector size. Each row corresponds to the k-dimensional label vector of a node. 
%
% fevo stores the evolution of graph energy at each iteration.
%
% Written by Amit Kumar K C, PhD at Wollman Lab (UCSD)
% Date: 5th May 2016
