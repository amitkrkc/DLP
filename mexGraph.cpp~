#include "mex.h"
#include "matrix.h"
#include "graph.h"

using namespace std;

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    //sanity checks
    if(nrhs!=5)
        mexErrMsgTxt("From mexGraph: invalid number of input arguments.\n");

    if(nlhs!=0)
        mexErrMsgTxt("From mexGraph: invalid number of outputs.\n");

    //associate parameter and perform sanity checks
    mxArray *params_in_m=mxDuplicateArray(prhs[0]);
    if(!mxIsDouble(params_in_m))
        mexErrMsgTxt("From mexGraph: invalid type of PARAMETERS. It must be a DOUBLE.\n");
    if(mxGetN(params_in_m)!=11)     // check the # of columns
        mexErrMsgTxt("From mexGraph: invalid number of PARAMETERS. It must be 1x11 vector.\n");
    double *params=mxGetPr(params_in_m);            // get pointer to params


    // associate edge matrix and perform sanity check
    mxArray *edge_in_m=mxDuplicateArray(prhs[1]);
    if(!mxIsDouble(edge_in_m))
        mexErrMsgTxt("From mexGraph: invalid type of EDGE_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(edge_in_m)!=3)
        mexErrMsgTxt("From mexGraph: invalid number of EDGE_IN_M. It must be nx3 matrix where n is the number of edges.\n");
    double *edgelist=mxGetPr(edge_in_m);


    // associate nuclei and perform sanity check
    mxArray *nuclei_in_m=mxDuplicateArray(prhs[2]);
    if(mxGetClassID(nuclei_in_m)!=mxINT32_CLASS)
        mexErrMsgTxt("From mexGraph: invalid type of NUCLEI_IN_M. It must be a INT32.\n");
    if(mxGetN(nuclei_in_m)!=1)
        mexErrMsgTxt("From mexGraph: invalid number of NUCLEI_IN_M. It must be nx1 matrix where n is the number of nuclei.\n");
    int *nuclei=(int *)mxGetData(nuclei_in_m);


    // associate score and perform sanity check
    mxArray *score_in_m=mxDuplicateArray(prhs[3]);
    if(!mxIsDouble(score_in_m))
        mexErrMsgTxt("From mexGraph: invalid type of SCORE_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(score_in_m)!=1)
        mexErrMsgTxt("From mexGraph: invalid number of SCORE_IN_M. It must be nx1 matrix where n is the number of superpixels.\n");
	double *score=mxGetPr(score_in_m);

    // associate lut and perform sanity check
    mxArray *lut_in_m=mxDuplicateArray(prhs[4]);
    if(!mxIsDouble(lut_in_m))
        mexErrMsgTxt("From mexGraph: invalid type of LUT_IN_M. It must be a DOUBLE.\n");
    if(mxGetN(lut_in_m)!=3)
        mexErrMsgTxt("From mexGraph: invalid number of LUT_IN_M. It must be nx3 matrix where n is the number of red cell superpixels.\n");
    int *lut=(int *)mxGetData(lut_in_m);



    // get the options
    options opts(params);

    //some more sanity checks
    if(opts.num_edges!=mxGetM(edge_in_m))
        mexErrMsgTxt("From mexGraph: number of edges in edgelist is and parameters do not match.\n");

    if(opts.lut_size!=mxGetM(lut_in_m))
        mexErrMsgTxt("From mexGraph: lut_size in parameters does not match with lut_in_m.\n");

     if(opts.num_nuclei!=mxGetM(nuclei_in_m))
        mexErrMsgTxt("From mexGraph: num_nuclei in parameters does not match with nuclei_in_m.\n");

     if(opts.num_nodes!=mxGetM(score_in_m))
        mexErrMsgTxt("From mexGraph: num_nodes in parameters does not match with score_in_m.\n");



    OUTPUT<<"BEFORE CALLING GRAPH"<<endl<<"Nuclei:";
    for(int i=0;i<opts.num_nuclei;i++)
        OUTPUT<<nuclei[i]<<"\t";
    OUTPUT<<endl;

    OUTPUT<<"Score:";
    for(int i=0;i<opts.num_nodes;i++)
        OUTPUT<<score[i]<<"\t";
    OUTPUT<<endl;


    ////////////////////////////////////////////////// ///////////////////////////////////////////
    Graph g;


    opts.display();
    g.add_node_and_edges(edgelist, nuclei, score, opts);
    g.display();


    /////////////////////////////// Gather labels
    vertex_descriptor p=0;
    OUTPUT<<"Gathering labels at node: "<<p<<endl;
    vector_double v_p(opts.label_vector_size+1);            //+1 for false positive
    double c_p;

    g.gather_labels_node(v_p, c_p, p, opts);

    for(int i=0;i<v_p.size();i++)
        OUTPUT<<v_p[i]<<"\t";
    cout<<endl;

    OUTPUT<<"Scheduling..."<<endl;
    vector<vertex_vector> J_all;
    g.schedule_nodes(J_all, 2);

    for(int i=0;i<J_all.size(); i++)
    {
        for(int j=0;j<J_all[i].size(); j++)
            cout<<J_all[i][j]<<"\t";
        cout<<endl;
    }
    cout<<endl;

    vector_double y_gt(opts.label_vector_size+1,0.0);
    vector_double z=g.solve_dc_nodewise(p,y_gt, opts);

    cout<<"Nodewise label:"<<endl;
    display<double>(&z[0], z.size());


    cout<<"Energy before:"<<g.compute_energy()<<endl;

    g.set_label(p, z);
    cout<<"Energy after:"<<g.compute_energy()<<endl;


    cout<<"argmax of label at p:"<<g.argmax(p)<<endl;
}
