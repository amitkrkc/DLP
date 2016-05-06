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








// check if the vertex is in the (nuclei) list
inline bool is_nuclei(const vector<bool>& nuclei_list, int q)
{
    return nuclei_list[q];
}


// add node(s) and edges
void add_node_and_edges(graph_t* g, double *edgelist, int *nuclei, double *score, double* params)
{

	struct timeval t_begin, t_end;
	gettimeofday(&t_begin,NULL);


//     params=[m k+1 1e-5 10 ~true, ~true 12]; % num_nodes, label_vector_size, tol, maxIter, verbose, use_parallel_processing, num_processors, num_edges
    const int num_nodes=(int)params[0];
    const int label_vector_size=(int)params[1];
    const int verbose=(int)params[4];
    const int num_edges=(int)params[7];
    const int num_nuclei=(int)params[8];

	if(verbose>1)
		OUTPUT<<"\n\t creating graph...";

    vertex_vector vlist(num_nodes);     //pre-allocating vlist (for speed?)

    // add all nuclei to the set
    // nuclei_list is a num_nodes x 1 array of booleans to store if a super-pixel is a nucleus or not.
    // it offers a constant time look-up later
    vector<bool> nuclei_list(num_nodes, false);
    for(int k=0;k<num_nuclei;k++)
    {
        int q=nuclei[k]-1;      //note the -1 because of C convention
        nuclei_list[q]=true;
    }

    	// add nodes to the graph
	for(int k=0;k<num_nodes;k++)
	{
        	// initialize the label vector
		vector_double label(label_vector_size,0.0);
		if(is_nuclei(nuclei_list, k))                //should be efficient than the previous version. Compared to earlier set based solution, this solution is O(1)
			label[k+1]=1.0;         //if it is a nucleus, assigns a unity vector, note k+1 instead of k. This is because of the fact that 0-th position is for the false positive
		else
			random_label(label);    // assigns a random label vector

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

	for(int k=0;k<num_edges;k++)
	{
		int s=(int)(*srclist++)-1;		//-1 because of c convention
		int t=(int)(*destlist++)-1;
		double w=*weightlist++;
		add_edge(vlist[s],vlist[t],arc(w),*g);
	}

    	gettimeofday(&t_end,NULL);
	double elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
	if(verbose>1)
		cout<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
}

/*
    solve label propagation at node p
    INPUTS:
        p: the node p at which label propagation is to be performed
        g: a (const) pointer to the graph
        Y: a (const) pointer to the label assignment matrix
        Y_gt: a const pointer to the label matrix obtained by using the red cells
        params: parameters
    OUTPUT:
        a vector double that corresponds to the updated label of the p-th node
*/
vector_double solve_dc_nodewise(vertex_descriptor p,const graph_t *g, const double *params)        //REMINDER: const double *Y means that data is const but pointer is not
{


    // REMINDER:     params=[num_nodes, label_vector_size, tol, maxIter, verbose, use_parallel_processing, num_processors, num_edges]
    const int label_vector_size=(int)params[1];
    const int verbose=(int)params[4];
    const double tol=params[2];

    if(verbose==2)
        cout<<"Inside solve_dc_nodewise function\n";


    // initialize v
    vector_double v=(*g)[p].label_gt;

    vector_double v(label_vector_size);
    memcpy(&v[0], Y+p, label_vector_size*sizeof(double));

    v[0]+=(*g)[p].score;

	// gather labels from the neighbors
	//out-edges of p
	out_edge_iterator oi,oi_end;
	double beta=(*g)[p].score;
	for(tie(oi,oi_end)=out_edges(p,*g);oi!=oi_end;oi++)
	{
		vertex_descriptor t=target(*oi,*g);
		vector_double zj=(*g)[t].label;
		double weight=(*g)[*oi].weight;
		if (verbose==2)     //2--> full verbose
		{
			cout<<"source:"<<p<<"\t target:"<<t<<"\t weight:"<<weight<<endl;
			display<double>(&zj[0], zj.size());
		}
		beta+=weight;
		for(uint i=0;i<label_vector_size;i++)
            		v[i]+=weight*zj[i];
	}

	if(verbose==2)
	{
		cout<<"-----------------------------------------\n";
		cout<<"beta:"<<beta<<endl;
		display<double>(&v[0], v.size());
		cout<<"-----------------------------------------\n";
	}

	// solve the problem
	vector_double z(label_vector_size);
	solve_quadratic_simplex(z, beta	, v, tol, verbose);	// TODO: check if beta and c are consistent

	if(verbose>1)
	{
		cout<<"Returning from optimization\t optimal label:";
		display<double>(&z[0], z.size());
	}
	return z;
}

void getLabels(graph_t *g, double* Y)
{
	int num_nodes=num_vertices(*g);
	int label_vector_size=(*g)[0].label.size();

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

	int num_nodes=num_vertices(*g);
	int label_vector_size=(*g)[0].label.size();

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
		//(*g)[i].label_gt=temp;
		ptr1++;
	}
}

// computes energy of the graph
double compute_energy(graph_t *g)
{
	double f=1.0;  //1 due to the term in \Delta

	vertex_iterator vi, vi_end;

	uint label_vector_size=(*g)[0].label.size();

    	// background label
	vector_double l0(label_vector_size, 0.0);
	l0[0]=1;


	for(tie(vi, vi_end)=vertices(*g); vi!=vi_end; vi++)
	{
		double fp=(*g)[*vi].score*norm2((*g)[*vi].label,l0);

		//subtract the ground-truth energy
		double gt_energy=0.0;
        for(uint i=0;i<label_vector_size;i++)
            gt_energy+=(*g)[*vi].label_gt[i]*(*g)[*vi].label[i];
        fp=fp-2*gt_energy;

		// out edges
		out_edge_iterator oi,oi_end;
		for(tie(oi,oi_end)=out_edges(*vi,*g);oi!=oi_end;oi++)
		{
			vertex_descriptor t=target(*oi,*g);
			fp+=(*g)[*oi].weight*norm2((*g)[*vi].label,(*g)[t].label);
		}
		f+=fp;
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
vector_double performLabelPropagation(graph_t *g, double *params, int *lut) //int maxIter=1000, double tol=1e-6, bool verbose=false, bool USE_PARALLEL_PROCESSING=true, int NUM_PROCESSORS=1)
{
    	int num_nodes=(int)params[0];
	int label_vector_size1=(int)params[1];
	double tol=params[2];
	int maxIter=(int)params[3];
	int verbose=(int)params[4];
	bool USE_PARALLEL_PROCESSING=(bool)params[5];
	int NUM_PROCESSORS=(int)params[6];
    	int num_edges=(int) params[7];
    	int lut_size=(int)params[8];

    if(verbose)
    {
        cout<<"\nInside performLabelPropagation..."<<endl
            <<"\tParameters:"<<endl
            <<"\t\t num_nodes:"<<num_nodes<<endl
            <<"\t\t label_vector_size1:"<<label_vector_size1<<endl
            <<"\t\t tol:"<<tol<<endl
            <<"\t\t maxIter:"<<maxIter<<endl
            <<"\t\t verbose:"<<verbose<<endl
            <<"\t\t USE_PARALLEL_PROCESSING:"<<USE_PARALLEL_PROCESSING<<endl
            <<"\t\t num_edges:"<<num_edges<<endl
            <<"\t\t lut_size:"<<lut_size<<endl;
    }
    //if(USE_PARALLEL_PROCESSING)
    //    omp_set_num_threads(NUM_PROCESSORS);
	vertex_iterator vi, vi_end;
	vector_double fevo;
	double max_change=-100000000000; //initialized to a very small number
	int iter=0;
	int label_vector_size=(*g)[0].label.size();

    	assert(label_vector_size==label_vector_size1);

    //allocate memory for gt_label
    double *gt_label;
    gt_label=new double[num_nodes*label_vector_size];

    double *Y;
    Y=new double[num_nodes*label_vector_size];

    float *params1;
    params1=new float[4];
    params1[0]=num_nodes;
    params1[1]=label_vector_size;
    params1[2]=lut_size;
    params1[3]=verbose;

	// schedule the vertices once
	vector<vertex_vector> J_all;
	schedule_nodes(J_all, g, NUM_PROCESSORS);

	// loop until convergence of the labels
	while(1)
	{
		iter++;

        	struct timeval t_begin, t_end;
        	gettimeofday(&t_begin,NULL);

        	if(verbose>1)
            		cout<<"\t computing gt_label...";

        	// compute gt-label
        	getLabels(g, Y);
        	compute_gt_label<double>(gt_label, Y, lut, params1);
        	if(verbose>1)
            		cout<<"done"<<endl;

        	if(verbose>1)
            		cout<<"\t Updating gt_label in the graph...";
        	update_gt_label_graph(g, gt_label);

        	if(verbose>1)
            		cout<<"done"<<endl;

        	// perform node-wise label propagation
		for(int ii=0;ii<J_all.size();ii++)
		{
			vertex_vector J=J_all[ii];

			// This loop CAN BE run in parallel
			for(int i=0;i<J.size();i++)
            		{
                		vertex_descriptor p=J[i];
                		vector_double zafter=solve_dc_nodewise(p,g,params);
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

	        if(verbose>0)
        		cout<<"iteration:"<<iter<<"\tgraph energy:"<<eval<<"\t max change:"<<max_change<<"\t fun change:"<<fun_change<<"\ttime taken:"<<elapsed_ms<<endl;

		if(iter>maxIter || fun_change<tol || max_change<tol)
			break;
	}

	// free the memory
	delete[] gt_label;
    delete[] Y;
    delete[] params1;

    if(verbose>0)
        cout<<endl;

	return fevo;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    struct timeval t_begin, t_end;              // to estimate time taken by each step
    double elapsed_ms;                          // time taken by each step

	mxArray *edge_in_m, /* *gt_in_m, */ *nuclei_in_m, *score_in_m, *lut_in_m;
	double *edge, *nuclei, *score;
	mxArray *params_in_m;
	double *params;
	int *lut;

	// associate inputs
	params_in_m=mxDuplicateArray(prhs[0]);
	edge_in_m=mxDuplicateArray(prhs[1]);
	//gt_in_m=mxDuplicateArray(prhs[2]);
	nuclei_in_m=mxDuplicateArray(prhs[2]);
	score_in_m=mxDuplicateArray(prhs[3]);
    lut_in_m=mxDuplicateArray(prhs[4]);

    //get pointers
    params=mxGetPr(params_in_m);
    edge=mxGetPr(edge_in_m);
    //gt=mxGetPr(gt_in_m);
	nuclei=mxGetPr(nuclei_in_m);
	score=mxGetPr(score_in_m);
    lut=(int *)mxGetData(lut_in_m);

	//extract parameters
	int num_nodes=(int)params[0];
	int label_vector_size=(int)params[1];
	double tol=params[2];
	int maxIter=(int)params[3];
	int verbose=(int)params[4];
	bool USE_PARALLEL_PROCESSING=(bool)params[5];
	int NUM_PROCESSORS=(int)params[6];
    int num_edges=(int) params[7];
	int lut_size=(int) params[8];
    if(verbose>0)
    {
        cout<<"parameters:\n\tnum_nodes:"<<num_nodes<<endl
            <<"\tlabel vector size:"<<label_vector_size<<endl
            <<"\ttol:"<<tol<<endl
            <<"\tmaxiter:"<<maxIter<<endl
            <<"\tverbose:"<<verbose<<endl
            <<"\tuse parallel:"<<USE_PARALLEL_PROCESSING<<endl
            <<"\tnum processors:"<<NUM_PROCESSORS<<endl
            <<"\tnum edges:"<<num_edges<<endl
	    <<"\tlut size:"<<lut_size<<endl;
    }
    // ensure that the edgelist have same number of edges
	assert(num_edges==(int) mxGetM(edge_in_m));


	//create graph
	if(verbose>0)
	{
        cout<<"Creating graph...";
	}

    gettimeofday(&t_begin,NULL);
	graph_t g;
    add_node_and_edges(&g, edge, nuclei, score, params); //add_node_and_edges(graph_t* g, double *edgelist, double *nuclei, double *score, double* params)

    gettimeofday(&t_end,NULL);
    elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
    if(verbose>0)
    {
        cout<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
    }

    if(verbose>0)
	{
        cout<<"Performing label propagation\n";
	}

    gettimeofday(&t_begin,NULL);
    vector_double fevo=performLabelPropagation(&g, params, lut);
    gettimeofday(&t_end,NULL);
    elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
    if(verbose>0)
    {
        cout<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
    }

	//retrieve labels from the graph
	if(verbose>0)
	{
        cout<<"Retrieving labels...";
	}
    gettimeofday(&t_begin,NULL);

	mxArray *Y_out_m;
	double *Y;
	Y_out_m=plhs[0]=mxCreateDoubleMatrix(num_nodes, label_vector_size, mxREAL);
	Y=mxGetPr(Y_out_m);
	getLabels(&g, Y);

    gettimeofday(&t_end,NULL);
    elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
    if(verbose>0)
    {
        cout<<"done...took "<<elapsed_ms<<" miliseconds"<<endl;
    }


	//retrieve energy evolution
	if(verbose>0)
	{
        cout<<"Retrieving energy evolution...";
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

