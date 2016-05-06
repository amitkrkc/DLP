#include "mex.h"
#include "matrix.h"
#include "dlp.h"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    //sanity checks
    if(nrhs!=5 || nlhs!=2)
        mexErrMsgTxt("From mexDLP: Incorrect format. Correct format is [Y,fevo]=mexDLP(params, edgelist, nuclei, score, lut) \n");

    //associate parameter and perform sanity checks
    mxArray *params_in_m=mxDuplicateArray(prhs[0]);
    if(!mxIsDouble(params_in_m))
        mexErrMsgTxt("From mexDLP: invalid type of PARAMETERS. It must be a DOUBLE.\n");
    if(mxGetN(params_in_m)!=11)     // check the # of columns
        mexErrMsgTxt("From mexDLP: invalid number of PARAMETERS. It must be 1x11 vector.\n");
    double *params=mxGetPr(params_in_m);            // get pointer to params


    // associate edge matrix and perform sanity check
    mxArray *edge_in_m=mxDuplicateArray(prhs[1]);
    if(!mxIsDouble(edge_in_m))
        mexErrMsgTxt("From mexDLP: invalid type of EDGE_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(edge_in_m)!=3)
        mexErrMsgTxt("From mexDLP: invalid number of EDGE_IN_M. It must be nx3 matrix where n is the number of edges.\n");
    double *edgelist=mxGetPr(edge_in_m);


    // associate nuclei and perform sanity check
    mxArray *nuclei_in_m=mxDuplicateArray(prhs[2]);
    if(mxGetClassID(nuclei_in_m)!=mxINT32_CLASS)
        mexErrMsgTxt("From mexDLP: invalid type of NUCLEI_IN_M. It must be a INT32.\n");
    if(mxGetN(nuclei_in_m)!=1)
        mexErrMsgTxt("From mexDLP: invalid number of NUCLEI_IN_M. It must be nx1 matrix where n is the number of nuclei.\n");
    int *nuclei=(int *)mxGetData(nuclei_in_m);


    // associate score and perform sanity check
    mxArray *score_in_m=mxDuplicateArray(prhs[3]);
    if(!mxIsDouble(score_in_m))
        mexErrMsgTxt("From mexDLP: invalid type of SCORE_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(score_in_m)!=1)
        mexErrMsgTxt("From mexDLP: invalid number of SCORE_IN_M. It must be nx1 matrix where n is the number of superpixels.\n");
	double *score=mxGetPr(score_in_m);

    // associate lut and perform sanity check
    mxArray *lut_in_m=mxDuplicateArray(prhs[4]);
    if(!mxIsDouble(lut_in_m))
        mexErrMsgTxt("From mexDLP: invalid type of LUT_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(lut_in_m)!=3)
        mexErrMsgTxt("From mexDLP: invalid number of LUT_IN_M. It must be nx3 matrix where n is the number of red cell superpixels.\n");
    int *lut=(int *)mxGetData(lut_in_m);



    // get the options
    options opts(params);

    //some more sanity checks
    if(opts.num_edges!=mxGetM(edge_in_m))
        mexErrMsgTxt("From mexDLP: number of edges in edgelist is and parameters do not match.\n");

    if(opts.lut_size!=mxGetM(lut_in_m))
        mexErrMsgTxt("From mexDLP: lut_size in parameters does not match with lut_in_m.\n");

     if(opts.num_nuclei!=mxGetM(nuclei_in_m))
        mexErrMsgTxt("From mexDLP: num_nuclei in parameters does not match with nuclei_in_m.\n");

     if(opts.num_nodes!=mxGetM(score_in_m))
        mexErrMsgTxt("From mexDLP: num_nodes in parameters does not match with score_in_m.\n");



    ////////////////////////////////////////////////// ///////
    // create a DLP object
    DLP dlp(edgelist, nuclei, score, lut, opts);

    // perform label propagation
    vector_double fevo;
    dlp.perform_label_propagation(fevo, opts);

    plhs[0]=mxCreateDoubleMatrix(opts.num_nodes, opts.label_vector_size+1, mxREAL);
    double *Y=mxGetPr(plhs[0]);
    dlp.get_labels(Y);


    plhs[1]=mxCreateDoubleMatrix(fevo.size(), 1, mxREAL);
    double *evo=mxGetPr(plhs[1]);
    memcpy(evo, &fevo[0],sizeof(double)*fevo.size());

}
