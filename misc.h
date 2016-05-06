#ifndef MISC_H
#define MISC_H

#include <iostream>     // for cout
#include <cstring>      // for memcpy, memset, etc
#include <cmath>            //for fabs()
#include <vector>       // for std::vector
#include <cstdlib>      // for drand48()
#include <cassert>      //for assert

#ifndef EPS
#define EPS 0.000001
#endif

#ifndef OUTPUT
#define OUTPUT std::cout
#endif

// some useful type definitions
typedef std::vector<double> vector_double;
typedef unsigned int uint;


// structure to hold all options
struct options
{
    int num_nodes;
	int label_vector_size;
	double tol;
	int maxIter;
	int verbose;
	bool use_parallel_processing;
	int num_processors;
    int num_edges;
    int lut_size;
    int num_nuclei;
    bool is_train_phase;

    // constructor
    options(const double *params)
    {
        num_nodes=(int)params[0];
        label_vector_size=(int)params[1];
        tol=params[2];
        maxIter=(int)params[3];
        verbose=(int)params[4];
        use_parallel_processing=(bool)params[5];
        num_processors=(int)params[6];
        num_edges=(int) params[7];
        lut_size=(int)params[8];
        num_nuclei=(int)params[9];
        is_train_phase=(bool) params[10];
    }

    // display the options
    void display() const
    {
        OUTPUT<<"OPTIONS:"<<std::endl
            <<"\tnum_nodes:"<<num_nodes<<std::endl
            <<"\tlabel vector size:"<<label_vector_size<<std::endl
            <<"\ttol:"<<tol<<std::endl
            <<"\tmaxiter:"<<maxIter<<std::endl
            <<"\tverbose:"<<verbose<<std::endl
            <<"\tuse parallel:"<<use_parallel_processing<<std::endl
            <<"\tnum processors:"<<num_processors<<std::endl
            <<"\tnum edges:"<<num_edges<<std::endl
            <<"\tlut size:"<<lut_size<<std::endl
            <<"\tnum_nuclei:"<<num_nuclei<<std::endl
            <<"\tis_train_phase:"<<is_train_phase<<std::endl;
    }
};

// check if the vertex is in the (nuclei) list
inline bool is_nuclei(const std::vector<bool>& nuclei_list, int q)
{
    return nuclei_list[q];
}


// computes argmax of the i-th row of Y
int argmax(double* Y, int i, int num_rows, int num_cols)
{
	assert(num_cols>=0);
	assert(num_rows>=0);
	assert(i>=0 && i<num_rows);      //sanity checks

    // initialize max_id
	int max_id=0;
	double *ptr=Y+i;
	double max_val=*ptr;

    // iterate through each element of the i-th row to find the maximum element
	for(int j=1;j<num_cols;j++)
	{
		ptr+=num_rows; //point to next item
		if(*ptr>max_val)
		{
			max_val=*ptr;
			max_id=j;
		}
	}
	return max_id;
}


//sets the i-th row to unit vector with 1 at the maxid index and zeros everywhere
// Y_gt is the output matrix of which i-th row is set to the unit vector
// val--> value to be set at the maxid-th column of the i-th row of Y_gt, it's default value is 1
void set_row(double *Y_gt, const int& i, const int& maxid, const int& num_rows, const int& num_cols, const double& val=1.0)
{
    assert(maxid>=0 && maxid<num_cols); // sanity checks
    assert(i>=0 && i<num_rows);


    double *ptr_out1=Y_gt+i;        //point to the i-th row of Y_gt

    //set the whole row to zero
    for(int j=0;j<num_cols;j++)
    {
        *ptr_out1=0.0;          //reset it to zero
        ptr_out1+=num_rows;     //advance the pointer to point to the next column
    }
    //set the maxid-th element to val (default=1)
    *(Y_gt+maxid*num_rows+i)=val;
}



// get the p-th row of Y_gt to y_gt
void get_row(vector_double& y_gt, int p, double *Y_gt, const options& opts)
{
    assert(p>=0 && p<opts.num_nodes); // sanity checks
    assert(y_gt.size()==opts.label_vector_size+1);      //+1 due to false positive label


    double *ptr=Y_gt+p;        //point to the p-th row of Y_gt

    //set the whole row to zero
    for(int j=0;j<y_gt.size();j++)
    {
        y_gt[j]=*ptr;           // store the value
        ptr+=opts.num_nodes;    //advance the pointer to point to the next column
    }
}




// computes square of an element
template<class T>
inline T SQUARE(const T& a)
{
	return a*a;
}


//computes L-2 norm of two vectors a and b
double norm2(const vector_double& a, const vector_double& b)
{
    assert(a.size()==b.size());

	double err=0.0;
	for(uint i=0;i<a.size();i++)
	{
		err+=SQUARE(a[i]-b[i]);
    }
	return err;
}

// computes dot product of two vectors
double dot_product(const vector_double& a, const vector_double& b)
{
    assert(a.size()==b.size());
    double sum=0.0;
    for(uint i=0;i<a.size();i++)
    {
        sum+=a[i]*b[i];
    }
    return sum;
}


// a random label vector of size m
void random_label(vector_double& x)
{
	double sum=EPS;

	for(uint i=0;i<x.size();i++)
	{
		x[i]=drand48();     //generate a random number [0,1)
		sum+=x[i];
	}
	for(uint i=0;i<x.size();i++)
	{
		x[i]/=sum;
	}
}

// display an array x
template<class T>
void display(T* x, int len)
{
	for(int i=0;i<len;i++)
	{
		OUTPUT<<x[i]<<"\t";
    }
	OUTPUT<<std::endl;
}

#endif
