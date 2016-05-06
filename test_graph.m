clear all;
close all;
clc;

num_nodes=10;
label_vector_size=5;
tol=1e-4;
maxIter=10;
verbose=2;
USE_PARALLEL_PROCESSING=0;
NUM_PROCESSORS=1;
lut_size=3;
is_train_phase=false;


A=triu(ones(num_nodes)-eye(num_nodes));
[rr,cc]=find(A>0);
edgelist=[rr cc randn(length(rr),1)];
num_edges=size(edgelist,1);

nuclei=int32([1;3;5]);
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
        num_nuclei,...
        is_train_phase
       ]

score=randn(num_nodes,1);

lut=randi(3,3);

% double* edgelist=mxGetPr(prhs[0]);
%    int* nuclei=(int *)mxGetData(prhs[1]);
%    double* score=mxGetPr(prhs[2]);
%    double* params=mxGetPr(prhs[3]);
%
%    Graph g;
%    options opts(params);
%    g.add_node_and_edges(edgelist, nuclei, score, opts);
%
mexGraph(params, edgelist, nuclei, score, lut );


disp('done');
