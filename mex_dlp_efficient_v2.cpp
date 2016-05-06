#include "misc.h"		// contains miscellaneous functions
#include "simplex.h"		// contains the function to solve the non-convex quadratic problem

#include <iostream>
#include <vector>
//#include <fstream>
#include <algorithm> 	//for copy
//#include <list>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <cstdlib>   //for drand48
#include <cmath>

//#include "omp.h"
#include <sys/time.h>
#include <cassert>
#include <cstring>

#include "mex.h"            //for mex
#include "matrix.h"


//#include "lutdata.cpp"

#include <set> //for set

#ifndef OUTPUT
#define OUTPUT std::cout
#endif
//#define USE_PARALLEL_PROCESSING 0
//#define NUM_PROCESSORS 10
#define EPS 1e-6

using namespace std;
using namespace boost;


//---------------------------------------------------
struct node
{
	vector_double label;
    vector_double label_gt;     //=-b_i*y_i^n/4/R+(1-b_i)*y_i^bar/4/Rbar
    double score;

    // constructors
	node(){}
	node(const vector_double& l, const double& s)
	{
		label=l;
		score=s;
		label_gt.clear();
	}
	node(const vector_double& l, const double& s, const vector_double& temp1)
	{
		label=l;
        label_gt=temp1;
        score=s;
	}
};

struct arc
{
	double weight;
	arc(double w=0.0)
	{
		weight=w;
	}
};
//------------------------------------------------------


//------------------typedef-------------------------
typedef adjacency_list<listS,vecS,undirectedS,node,arc > graph_t;
typedef graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
typedef graph_traits<graph_t>::vertex_iterator vertex_iterator;
typedef graph_traits<graph_t>::edge_iterator edge_iterator;
typedef graph_traits<graph_t>::in_edge_iterator in_edge_iterator;
typedef graph_traits<graph_t>::out_edge_iterator out_edge_iterator;
typedef graph_traits<graph_t>::edge_descriptor edge_descriptor;
typedef vector<vertex_descriptor> vertex_vector;
typedef list<vertex_descriptor> vertex_list;
typedef pair<int,int> edge_t;
typedef set<int> set_int;
typedef set<vertex_descriptor> vertex_set;


// add node(s) and edges
void add_node_and_edges(graph_t* g, double *edgelist, int *nuclei, double *score, const options& opts)
{
	struct timeval t_begin, t_end;
	gettimeofday(&t_begin,NULL);


//     params=[m k+1 1e-5 10 ~true, ~true 12]; % num_nodes, label_vector_size, tol, maxIter, verbose, use_parallel_processing, num_processors, num_edges
    const int num_nodes=(int)params[0];
    const int label_vector_size=(int)params[1];
    const int verbose=(int)params[4];
    const int num_edges=(int)params[7];
    const int num_nuclei=(int)params[8];

	if(opts.verbose>1)
		OUTPUT<<"\n\t creating graph...";

    vertex_vector vlist(opts.num_nodes);     //pre-allocating vlist (for speed?)

    // add all nuclei to the set
    // nuclei_list is a num_nodes x 1 array of booleans to store if a super-pixel is a nucleus or not.
    // it offers a constant time look-up later
    vector<bool> nuclei_list(opts.num_nodes, false);
    for(int k=0;k<opts.num_nuclei;k++)
    {
        int q=nuclei[k]-1;      //note the -1 because of C convention
        nuclei_list[q]=true;
    }

    	// add nodes to the graph
	for(int k=0;k<opts.num_nodes;k++)
	{
        	// initialize the label vector
		vector_double label(opts.label_vector_size,0.0);
		if(is_nuclei(nuclei_list, k))       //should be efficient than the previous version. Compared to earlier set based solution, this solution is O(1). is_nuclei() is defined in misc.h
			label[k+1]=1.0;         //if it is a nucleus, assigns a unity vector, note k+1 instead of k. This is because of the fact that 0-th position is for the false positive
		else
			random_label(label);    // assigns a random label vector, defined in misc.h

		// score s=alpha*p_k-beta*(1-p_k)
        double s=score[k];

        // add vertex
		node n(label,s);
		vertex_descriptor v=add_vertex(n,*g);
		vlist[k]=v;                             // replaced the old version: vlist.push_back(v), TODO: check if vlist[k]=v is faster than vlist.push_back(v)
	}

    //add edges
    double *srclist=edgelist;
    double *destlist=edgelist+num_edges;
    double *weightlist=edgelist+2*num_edges;

	for(int k=0;k<opts.num_edges;k++)
	{
		int s=(int)(*srclist++)-1;		//-1 because of c convention
		int t=(int)(*destlist++)-1;
		double w=*weightlist++;
		add_edge(vlist[s],vlist[t],arc(w),*g);
	}

    gettimeofday(&t_end,NULL);
	double elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
	if(opts.verbose>1)
		OUTPUT<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
}


// solve label propagation at node p
vector_double solve_dc_nodewise(vertex_descriptor p,const graph_t *g, const options& opts)       // REMINDER: const type *x: data is constant but pointer is not
{
    if(opts.verbose==2)
        OUTPUT<<"Inside solve_dc_nodewise function\n";

    // v_p=y_p^red+s_p*l_0+\sum_j w_pj^eff*y_j
    // initialize v
    vector_double v_p=(*g)[p].label_gt;   //y_p^red
    v_p[0]+=(*g)[p].score;                 // add s_p*l_0, since l_0=[1,0,0,...,0], only first index is added

	// gather labels from the neighbors

	// initialize c_p
	double c_p=(*g)[p].score;

	//edges of p
	edge_iterator ei,ei_end;
	for(tie(ei,ei_end)=edges(p,*g);ei!=ei_end;ei++)
	{
		vertex_descriptor t=target(*ei,*g);
		vector_double y_j=(*g)[t].label;
		double w_pj=(*g)[*ei].weight;

		if (opts.verbose==2)     //2--> full verbose
		{
			OUTPUT<<"source:"<<p<<"\t target:"<<t<<"\t weight:"<<w_pj<<endl;
			display<double>(&y_j[0], y_j.size());
		}

		c_p+=w_pj;      // add w_pj

		// gather labels
		for(uint i=0;i<opts.label_vector_size;i++)
            v_p[i]+=w_pj*y_j[i];
	}

	if(opts.verbose==2)
	{
		OUTPUT<<"-----------------------------------------\n";
		OUTPUT<<"c_p:"<<c_p<<endl;
		display<double>(&v[0], v.size());
		OUTPUT<<"-----------------------------------------\n";
	}

	// solve the quadratic problem at node p
	vector_double z(label_vector_size);
	solve_quadratic_simplex(z, c_p	, v_p, opts);	// z contains the updated label of node p

	if(opts.verbose>1)
	{
		OUTPUT<<"Returning from optimization\t optimal label:";
		display<double>(&z[0], z.size());
	}
	return z;
}

void getLabels(graph_t *g, double* Y)
{
	const int num_nodes=num_vertices(*g);
	const int label_vector_size=(*g)[0].label.size();

	double *ptr1, *ptr2;
	ptr1=Y;
	for(int i=0;i<num_nodes;i++)
	{
		vector_double temp=(*g)[i].label;
		ptr2=ptr1;
		for(int j=0;j<label_vector_size;j++)
		{
			*ptr2=temp[j];
			ptr2+=num_nodes;
		}
		ptr1++;
	}
}

void update_gt_label_graph(graph_t *g, double* gt_label)
{

	const int num_nodes=num_vertices(*g);
	const int label_vector_size=(*g)[0].label.size();

	double *ptr1, *ptr2;
	ptr1=gt_label;
	for(int i=0;i<num_nodes;i++)
	{
		ptr2=ptr1;
		//vector_double temp(label_vector_size,0.0);
		for(int j=0;j<label_vector_size;j++)
		{
			(*g)[i].label_gt[j]=*ptr2;
			ptr2+=num_nodes;
		}
		ptr1++;
	}
}

// computes energy of the graph
double compute_energy(graph_t *g)
{
	// total energy f, initialized to zero
	double f=0.0;

    // label vector size
    const uint label_vector_size=(*g)[0].label.size();

    // false positive label
	vector_double l0(label_vector_size, 0.0);
	l0[0]=1.0;

    // iterate through all vertices
    vertex_iterator vi, vi_end;
	for(tie(vi, vi_end)=vertices(*g); vi!=vi_end; vi++)
	{
        vertex_descriptor p=*vi;

        // label of the p-th vertex
        vector_double y_p=(*g)[p].label;

        // dot product of y_p^red and y_p
        double term1=dot_product((*g)[p].label_gt, y_p);

        // L2 norm of y_p-L_0
        double term2=0.5*(g*)[p].score*norm2(y_p, l_0);

        // label dissimilarity between the edges
        double term3=0.0;
        edge_iterator ei, ei_end;
        for(tie(ei, ei_end)=edges(*g, p); ei!=ei_end; ei++)
        {
            vertex_descriptor j=target(*ei,*g);     //node j
            vector_double y_j=(*g)[j].label;        // label of the j-th node
            double w_pj=(*g)[*ei].weight;           // weight of the (p,j) edge
            term3+=w_pj*norm2(y_p, y_j);            // label discrepancy between y_p and y_j
        }
        term3*=0.5;

        f+=term1+term2+term3;
	}
	return f;
}

// computes argmax of the label distribution
uint argmax(const vector_double& y)			//argmax(double* Y, int i, int num_rows, int num_cols)
{
//    return arg_max<double>(y);
    uint maxid=0;
    double maxval=y[0];
    for(uint i=1;i<y.size(); i++)
    {
        if(y[i]>maxval)
        {
            maxval=y[i];
            maxid=i;
        }
    }
    return maxid;
}

// performs label propagation at each edge until convergence
//return the graph energy evolution
vector_double performLabelPropagation(graph_t *g, int *lut, const options& opts)
{
    if(verbose)
    {
        OUTPUT<<"\nInside performLabelPropagation..."<<endl;
        opts.display();
    }

	vertex_iterator vi, vi_end;
	vector_double fevo;
	const double max_change=-100000000000; //initialized to a very small number
	int iter=0;

    // sanity check
    assert((*g)[0].label.size()==opts.label_vector_size);

    //allocate memory for gt_label
    double *gt_label;
    gt_label=new double[opts.num_nodes*opts.label_vector_size];

    // allocate memory for label assignment matrix Y
    double *Y;
    Y=new double[opts.num_nodes*opts.label_vector_size];


	// schedule the vertices once
	vector<vertex_vector> J_all;
	schedule_nodes(J_all, g, opts.NUM_PROCESSORS);

	// loop until convergence of the labels
	while(1)
	{
		iter++;

        struct timeval t_begin, t_end;
        gettimeofday(&t_begin,NULL);

        if(opts.verbose>1)
            OUTPUT<<"\t computing gt_label...";

        //---------------------------------TODO: check if this block of the code is the bottle-neck------------------------
        // compute gt-label
        getLabels(g, Y);
        compute_gt_label<double>(gt_label, Y, lut, opts);

        if(opts.verbose>1)
            OUTPUT<<"done"<<endl;

        if(opts.verbose>1)
            OUTPUT<<"\t Updating gt_label in the graph...";
        update_gt_label_graph(g, gt_label);

        if(opts.verbose>1)
            OUTPUT<<"done"<<endl;
        //---------------------------------TODO: check if this block of the code is the bottle-neck------------------------

        // perform node-wise label propagation
		for(int ii=0;ii<J_all.size();ii++)
		{
			vertex_vector J=J_all[ii];

			// This loop CAN BE run in parallel
			for(int i=0;i<J.size();i++)
            {
                vertex_descriptor p=J[i];
                vector_double zafter=solve_dc_nodewise(p,g,opts);
                max_change=MAX<double>(max_change, norm2((*g)[p].label, zafter));
                (*g)[p].label=zafter;
            }
		}

		double eval=compute_energy(g);

		double fun_change=100000000;
		if(fevo.size()>1)
	        	fun_change=fabs(eval-*fevo.rbegin());
		fevo.push_back(eval);

        gettimeofday(&t_end,NULL);
        double elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;

        if(opts.verbose>0)
            OUTPUT<<"iteration:"<<iter<<"\tgraph energy:"<<eval<<"\t max change:"<<max_change<<"\t fun change:"<<fun_change<<"\ttime taken:"<<elapsed_ms<<endl;

        // termination criteria
		if(iter>maxIter || fun_change<tol || max_change<tol)
			break;
	}

	// free the memory
	delete[] gt_label;
    delete[] Y;
    delete[] params1;

    if(opts.verbose>0)
        OUTPUT<<endl;

	return fevo;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //sanity checks
    if(nrhs!=5)
        mexErrorMsgTxt("From mexDLP: invalid number of input arguments.\n");

    if(nlhs!=1)
        mexErrorMsgTxt("From mexDLP: invalid number of outputs.\n");

    //associate parameter and perform sanity checks
    mxArray *params_in_m=mxDuplicateArray(prhs[0]);
    if(!mxIsDouble(params_in_m))
        mexErrorMsgTxt("From mexDLP: invalid type of PARAMETERS. It must be a DOUBLE.\n");
    if(mxGetN(params_in_m)!=10)
        mexErrorMsgTxt("From mexDLP: invalid number of PARAMETERS. It must be 1x10 vector.\n");
    double *params=mxGetPr(params_in_m);            // get pointer to params


    // associate edge matrix and perform sanity check
    mxArray *edge_in_m=mxDuplicateArray(prhs[1]);
    if(!mxIsDouble(edge_in_m))
        mexErrorMsgTxt("From mexDLP: invalid type of EDGE_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(edge_in_m)!=3)
        mexErrorMsgTxt("From mexDLP: invalid number of EDGE_IN_M. It must be nx3 matrix where n is the number of edges.\n");
    double *edge=mxGetPr(edge_in_m);


    // associate nuclei and perform sanity check
    mxArray *nuclei_in_m=mxDuplicateArray(prhs[2]);
    if(!mxIsDouble(nuclei_in_m))
        mexErrorMsgTxt("From mexDLP: invalid type of NUCLEI_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(nuclei_in_m)!=1)
        mexErrorMsgTxt("From mexDLP: invalid number of NUCLEI_IN_M. It must be nx1 matrix where n is the number of nuclei.\n");
    double *nuclei=mxGetPr(nuclei_in_m);


    // associate score and perform sanity check
    mxArray *score_in_m=mxDuplicateArray(prhs[3]);
    if(!mxIsDouble(score_in_m))
        mexErrorMsgTxt("From mexDLP: invalid type of SCORE_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(score_in_m)!=1)
        mexErrorMsgTxt("From mexDLP: invalid number of SCORE_IN_M. It must be nx1 matrix where n is the number of superpixels.\n");
	double *score=mxGetPr(score_in_m);

    // associate lut and perform sanity check
    mxArray *lut_in_m=mxDuplicateArray(prhs[4]);
    if(!mxIsDouble(lut_in_m))
        mexErrorMsgTxt("From mexDLP: invalid type of LUT_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(lut_in_m)!=3)
        mexErrorMsgTxt("From mexDLP: invalid number of LUT_IN_M. It must be nx3 matrix where n is the number of red cell superpixels.\n");
    double *lut=(int *)mxGetData(lut_in_m);


    struct timeval t_begin, t_end;              // to estimate time taken by each step
    double elapsed_ms;                          // time taken by each step

	//extract parameters
	const options opts(params);

    if(verbose>0)
    {
        opts.display();
    }

    // ensure that the edgelist have same number of edges
    if(opts.num_edges!=mxGetM(edge_in_m));
        mexErrorMsgTxt("From mexDLP: number of edges in PARAMETERS and in EDGE_IN_M do not match.\n");


	//create graph
	if(opts.verbose>0)
	{
        OUTPUT<<"Creating graph...";
	}

    gettimeofday(&t_begin,NULL);
	graph_t g;
    add_node_and_edges(&g, edge, nuclei, score, opts); //add_node_and_edges(graph_t* g, double *edgelist, double *nuclei, double *score, double* params)

    gettimeofday(&t_end,NULL);
    elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
    if(opts.verbose>0)
    {
        OUTPUT<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
    }

    if(opts.verbose>0)
	{
        OUTPUT<<"Performing label propagation\n";
	}

    gettimeofday(&t_begin,NULL);
    vector_double fevo=performLabelPropagation(&g, lut, opts);
    gettimeofday(&t_end,NULL);
    elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
    if(opts.verbose>0)
    {
        OUTPUT<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
    }

	//retrieve labels from the graph
	if(opts.verbose>0)
	{
        OUTPUT<<"Retrieving labels...";
	}
    gettimeofday(&t_begin,NULL);

	mxArray *Y_out_m;
	double *Y;
	Y_out_m=plhs[0]=mxCreateDoubleMatrix(num_nodes, label_vector_size, mxREAL);
	Y=mxGetPr(Y_out_m);
	getLabels(&g, Y);

    gettimeofday(&t_end,NULL);
    elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
    if(opts.verbose>0)
    {
        OUTPUT<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
    }


	//retrieve energy evolution
	if(opts.verbose>0)
	{
        OUTPUT<<"Retrieving energy evolution...";
	}

	gettimeofday(&t_begin,NULL);
	mxArray *evo_out_m;
	double *evo;
	evo_out_m=plhs[1]=mxCreateDoubleMatrix(fevo.size(), 1, mxREAL);
	evo=mxGetPr(evo_out_m);
	copy(fevo.begin(), fevo.end(), evo);
    getLabels(&g, Y);

    gettimeofday(&t_end,NULL);
    elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
    if(verbose>0)
    {
        cout<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
    }

	return;
}

