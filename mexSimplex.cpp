#include "mex.h"
#include "matrix.h"
#include "simplex.h"
// #include "misc.h"

using namespace std;


void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    double* in=mxGetPr(prhs[0]);
    const int vec_size=mxGetM(prhs[0]);
    
    vector_double v(vec_size);
    memcpy(&v[0], in, sizeof(double)*vec_size);
    
    double c=mxGetScalar(prhs[1]);
    
    simplex s(v, c);
    
    double* params=mxGetPr(prhs[2]);
    
    
    options opts(params);
    vector_double v_out=s.solve_quadratic_simplex(opts);
    
    
    plhs[0]=mxCreateDoubleMatrix(vec_size, 1, mxREAL);
    double* out=mxGetPr(plhs[0]);
    
    memcpy(out, &v_out[0], sizeof(double)*vec_size);
}
