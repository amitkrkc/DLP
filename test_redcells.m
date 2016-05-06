% clear all;close all;clc;




lut_size=5;

num_nodes=10;
label_vector_size=5;
tol=1e-4;
maxIter=10;
verbose=2;
USE_PARALLEL_PROCESSING=0;
NUM_PROCESSORS=1;

sp=int32([(1:2:num_nodes)' randi([1 2],lut_size,1) randi([1 100],lut_size,1)]);



A=triu(ones(num_nodes)-eye(num_nodes));
[rr,cc]=find(A>0);
edgelist=[rr cc randn(length(rr),1)];
num_edges=size(edgelist,1);

nuclei=[1;3;5];
num_nuclei=length(nuclei);

params=[
        num_nodes, ...
        label_vector_size,... 
        tol, ...
        maxIter,... 
        verbose, ...
        USE_PARALLEL_PROCESSING, ...
        NUM_PROCESSORS, ...
        num_edges, ...
        lut_size,...
        num_nuclei
       ]
   
 mexRedCells(params, sp);