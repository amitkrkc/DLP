#include "mex.h"            //for mex
#include "matrix.h"

#include "redcell.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if(nrhs!=2 || nlhs!=0)
    {
        mexErrMsgTxt("Invalid number of arguments. Correct syntax: mexRedCells(params, sp_in)");
    }

     //associate parameter and perform sanity checks
    mxArray *params_in_m=mxDuplicateArray(prhs[0]);
    if(!mxIsDouble(params_in_m))
        mexErrMsgTxt("From mexRedCells: invalid type of PARAMETERS. It must be a DOUBLE.\n");
    if(mxGetN(params_in_m)!=11)     // check the # of columns
        mexErrMsgTxt("From mexGraph: invalid number of PARAMETERS. It must be 1x11 vector.\n");
    double *params=mxGetPr(params_in_m);

    mxArray *sp_in_m=mxDuplicateArray(prhs[1]);
    if(mxGetClassID(sp_in_m)!=mxINT32_CLASS)
        mexErrMsgTxt("From mexRedCells: invalid type of SP_IN_M. It must be a INT32.\n");
    int *lut_data=(int *)mxGetData(sp_in_m);


    // get the options
    options opts(params);


    Redcell redcells(lut_data, opts);


    redcells.display();


//     cout<<"From mex function"<<endl
//         <<"num red cells:"<<redcells.num_red_cells<<endl
//         <<"num not red cells:"<<redcells.num_not_red_cells<<endl;

//    vector<bool> is_red=redcells.get_is_red();
//    display<bool>(&is_red[0], is_red.size());
}
