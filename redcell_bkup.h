#ifndef REDCELL_H
#define REDCELL_H

#include <algorithm>    //for sort
#include <cstring>      //for memcpy
#include <vector>
#include <iostream>
#include <set>
#include <cassert>
#include <map>
#include "mex.h"
#ifndef OUTPUT
#define OUTPUT std::cout
#endif

#include "misc.h"       //contains argmax and set_row functions

const int VERBOSE=0;

using namespace std;

//structure to hold the data of a red cell
// each red cell (that has a unique cell label), contains the total area, a set of superpixels
// and
struct red_cell
{

	int total_area;         //total area of the super pixels
	set<int> superpixels;   // list of superpixels

	map<int, int> maxid_area;	//first int-->max id, second int-->total area of
                                //the pixels that have the same max_id
};

typedef map<int, red_cell> REDCELL;




// used for sorting the red cells in decreasing order of total_area
struct greater_area_key
{
	bool operator() (const red_cell& a, const red_cell& b) const
	{
		return a.total_area>b.total_area;
	}
};

// insert the data (sp, curr_label, pix_area, max_id) to red cell
// (sp, curr_label, pix_area) corresponds to a row of input data
// max_id is computed as the argmax of the sp-th row of Y matrix
void insertToLUT(REDCELL& lut, const int& sp, const int& curr_label, const int& pix_area, const int& max_id)
{
    REDCELL::iterator it=lut.find(curr_label);

    //if the curr_label does not exist in the map, create a new one
    if(it==lut.end())
    {

        //create a new red_cell object temp
        red_cell temp;
        temp.total_area=pix_area;
        temp.superpixels.insert(sp);
        //temp.max_id.insert(max_id);
        temp.maxid_area[max_id]=pix_area;

        lut[curr_label]=temp;
    }
    //otherwise update the entry in the map
    else
    {
        it->second.total_area+=pix_area;    //add the total pixel area;
        it->second.superpixels.insert(sp);  //insert the superpixel_area

        //check if the max-id already exists
        map<int, int>::iterator it1=it->second.maxid_area.find(max_id);
        if(it1==it->second.maxid_area.end())
        {
            //if not found, create a new entry in the maxid_area map
            it->second.maxid_area[max_id]=pix_area;

        }
        else
        {
            //if found, update the pix_area
            it1->second+=pix_area;
        }
    }
}


// displays the content of each red cell
void display(REDCELL lut)
{
    for(REDCELL::iterator it=lut.begin(); it!=lut.end(); it++)
    {

        cout<<"============================================="<<endl
            <<"red cell label:"<<it->first<<"\t total area:"<<it->second.total_area
            <<"\tsuperpixels: {";

        //display all superpixels inside the red cell
        for(set<int>::iterator it1=it->second.superpixels.begin(); it1!=it->second.superpixels.end(); it1++)
            cout<<*it1<<"  ";
        cout<<"}"<<endl;

        // display max-id and area
        for(map<int, int>::iterator it1=it->second.maxid_area.begin(); it1!=it->second.maxid_area.end(); it1++)
        {
            cout<<"\t max id:"<<it1->first<<"\t tot area:"<<it1->second<<endl;
        }
        cout<<endl;
    }
}

// reads from the inputs lut_data and inserts into REDCELL
// lut_data is a Nx3 integer matrix where each row corresponds to (sp, curr_label, pix_area)
// lut is the output REDCELL structure
// set_not_red_cells: set of other superpixels (i.e., the superpixels that do not belong to the red cells)
void find_unique_labels(REDCELL& lut, set<int>& set_not_red_cells, int *lut_data, double *Y, int lut_size, int num_nodes, int label_vector_size)
{
    if(VERBOSE>0)
        OUTPUT<<"Inside find_unique_labels"<<endl;

    // first all all num_nodes in set_not_red_cells and initialize set_red_cells to empty
    // later, we will remove an element from set_not_red_cells if the elment belongs to the red cell and add it to the set_red_cells

//    set_red_cells.clear();            //DISABLED because set of superpixels belonging to the red cells is already known in REDCELL lut
    for(int i=0;i<num_nodes;i++)
        set_not_red_cells.insert(i);

    for(int i=0;i<lut_size; i++)
	{
		// read the current entry from the lut_data
		int sp=lut_data[i]-1;                       //NOTE: 1 is subtracted for C convention
		int curr_label=lut_data[i+lut_size];
		int pix_area=lut_data[i+2*lut_size];

		// remove sp from set_not_red_cells
		set_not_red_cells.erase(sp);

        //find the arg-max of the sp-th row of the label matrix
        int max_id=argmax(Y, sp, num_nodes, label_vector_size);


        if(VERBOSE>1)
        {
            OUTPUT<<"\n*************************************\n"
                  <<"sp:"<<sp<<"\t curr label:"<<curr_label<<"\t area:"<<pix_area<<"\t max id:"<<max_id<<endl;
        }

        insertToLUT(lut, sp, curr_label, pix_area, max_id);

        if(VERBOSE>1)
        {
            display(lut);
        }
    }
}


void compute_gt_label(double *Y_gt, double* Y, REDCELL lut, const set<int>& set_not_red_cells, const int& num_nodes, const int& label_vector_size)
{

    //pre-computing some factors
    const int num_not_red_cells=set_not_red_cells.size();
    const int num_red_cells=num_nodes-num_not_red_cells;

    const double factor_red_cells=-0.25/num_red_cells;
    const double factor_not_red_cells=1.0/num_red_cells/num_not_red_cells;

    // since map does not sort the red cells based on their area,
    // we first convert it to a vector and then invoke the sort function


    vector<red_cell> mapcopy;
    for(REDCELL::iterator it=lut.begin(); it!=lut.end(); it++)
    {
        mapcopy.push_back(it->second); //make_pair(it->first, it->second));
    }

    // sort the mapcopy based on the area
    std::sort(mapcopy.begin(), mapcopy.end(), greater_area_key() );

    if(VERBOSE>1)
        OUTPUT<<"AFTER SORTING\n";

    // copy Y_gt to zero initially
    // later, each row of Y_gt will be modified depending on whether the super-pixel belongs to the red-cell or not
    // we use memset because we will modify only few columns of the superpixels that belong to the red-cell
    memset(Y_gt, 0.0, sizeof(double)*num_nodes*label_vector_size);


    vector<int> max_labels; // required later for not-red-cell superpixels

    // for each red cell (which is sorted based on the area)
    for(vector<red_cell>::iterator it=mapcopy.begin(); it!=mapcopy.end(); it++)
    {
        if(VERBOSE>1)
        {
            OUTPUT<<it->total_area<<"\t superpixels:{";
            for(set<int>::iterator it1=it->superpixels.begin(); it1!=it->superpixels.end(); it1++)
                OUTPUT<<" "<<*it1;
            OUTPUT<<"}"<<endl;

        }

        //find the maxid that has the largest area
        int maxid=-1;
        int maxarea=-1;
        for(map<int, int>::iterator it1=it->maxid_area.begin(); it1!=it->maxid_area.end(); it1++)
        {
            if(maxarea<it1->second) //if the current max is smaller
            {
                maxarea=it1->second;
                maxid=it1->first;
            }
        }
        if(VERBOSE>1)
            OUTPUT<<"\tmaxid:"<<maxid<<"\t maxarea:"<<maxarea<<endl;

        // if maxid is not -1
        if(maxid!=-1)           //sanity check
        {
            // set all the red pixels to a vector with -1/4R at the maxid-th index and zero everywhere where R is the number of superpixels in the red cells
            for(set<int>::iterator it1=it->superpixels.begin(); it1!=it->superpixels.end(); it1++)
            {
                if(VERBOSE>1)
                    OUTPUT<<"\t setting "<<*it1<<" row to e-vector at "<<maxid<<endl;
                //set_row(Y_gt, *it1, maxid, num_nodes, label_vector_size);
                //DISABLED: set_row because Y_gt has already been memset to 0 before the beginning of the for loop
                // set_row also sets the row to zero. save extra loop,
                int sp=*it1;
                *(Y_gt+maxid*num_nodes+sp)=factor_red_cells;
                max_labels.push_back(maxid);
            }
        }
    }

    // for the superpixels that are not in the red-cells
    for(set<int>::iterator it=set_not_red_cells.begin(); it!=set_not_red_cells.end(); it++)
    {
        int not_sp=*it;
        for(uint j=0;j<max_labels.size();j++)
        {
            int lbl=max_labels[j];
            Y_gt[not_sp+num_nodes*lbl]=factor_not_red_cells;
        }
    }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int lut_size=(int)mxGetM(prhs[0]);
	int lut_dim=(int) mxGetN(prhs[0]);

	if(!mxIsInt32(prhs[0]))
		mexErrMsgTxt("Invalid data type. Input data type must be int.\n");

	if(lut_dim!=3)
		mexErrMsgTxt("Invalid size. Input matrix  must be nx3\n");

	int label_vector_size=(int)mxGetN(prhs[1]);
	int num_nodes=	(int)mxGetM(prhs[1]);

	if(!mxIsDouble(prhs[1]))
		mexErrMsgTxt("Invalid data type. Y must be double.\n");

    if(VERBOSE>0)
    {
        OUTPUT<<"lut_size:"<<lut_size<<"\t lut_dim:"<<lut_dim<<endl;
        OUTPUT<<"num nodes:"<<num_nodes<<"\tlabel vector size:"<<label_vector_size<<endl;
    }
	int *lut_data=(int *)mxGetData(prhs[0]);
	double *Y=mxGetPr(prhs[1]);

    //create output matrix and associate a pointer
    plhs[0]=mxCreateDoubleMatrix(num_nodes, label_vector_size, mxREAL);
    double *Y_gt=mxGetPr(plhs[0]);

    // read from the inputs and create REDCELL map
	REDCELL lut;
	set<int> set_not_red_cells;  // to store the set of superpixels that do not belong to the red cells

//	find_unique_labels(lut, lut_data, Y, lut_size, num_nodes, label_vector_size);
    find_unique_labels(lut, set_not_red_cells, lut_data, Y, lut_size, num_nodes, label_vector_size);


    if(VERBOSE>0)
    {
        OUTPUT<<"Red cells\n";
        display(lut);
    }



    compute_gt_label(Y_gt, Y, lut, set_not_red_cells, num_nodes, label_vector_size);

//    vector< red_cell > mapcopy;
//    for(REDCELL::iterator it=lut.begin(); it!=lut.end(); it++)
//    {
//        mapcopy.push_back(it->second);//make_pair(it->first, it->second));
//    }
//
//
//    std::sort(mapcopy.begin(), mapcopy.end(), greater_area_key() );
//
//    if(VERBOSE>1)
//        OUTPUT<<"AFTER SORTING\n";
//
//    plhs[0]=mxCreateDoubleMatrix(num_nodes, label_vector_size, mxREAL);
//    double *Y_gt=mxGetPr(plhs[0]);
//
//    // copy the Y matrix to Y_gt
//    memcpy(Y_gt, Y, sizeof(double)*num_nodes*label_vector_size);
//
//    for(vector<red_cell>::iterator it=mapcopy.begin(); it!=mapcopy.end(); it++)
//    {
//        if(VERBOSE>1)
//        {
//            OUTPUT<<it->total_area<<"\t superpixels:{";
//            for(set<int>::iterator it1=it->superpixels.begin(); it1!=it->superpixels.end(); it1++)
//                OUTPUT<<" "<<*it1;
//            OUTPUT<<"}"<<endl;
//
//        }
//
//        //find the maxid that has the largest area
//        int maxid=-1;
//        int maxarea=-1;
//        for(map<int, int>::iterator it1=it->maxid_area.begin(); it1!=it->maxid_area.end(); it1++)
//        {
//            if(maxarea<it1->second) //if the current max is smaller
//            {
//                maxarea=it1->second;
//                maxid=it1->first;
//            }
//        }
//        if(VERBOSE>1)
//            OUTPUT<<"\tmaxid:"<<maxid<<"\t maxarea:"<<maxarea<<endl;
//
//        if(maxid!=-1)           //sanity check
//        {
//            for(set<int>::iterator it1=it->superpixels.begin(); it1!=it->superpixels.end(); it1++)
//            {
//                if(VERBOSE>1)
//                    OUTPUT<<"\t setting "<<*it1<<" row to e-vector at "<<maxid<<endl;
//                set_row(Y_gt, *it1, maxid, num_nodes, label_vector_size);
//            }
//        }
//
////        std::sort(temp_vec.begin(), temp_vec.end(), greater_area_maxid() );
////
////        int maxid=*(temp_vec.begin());
////
////        for(set<int>::iterator it1=it->superpixels.begin(); it1!=it->superpixels.end(); it1++)
////        {
////
////        }
////
////
////        for(vector< pair<int, int> > ::iterator it1=temp_vec.begin(); it1!=temp_vec.end(); it1++)
////        {
////            cout<<"\t"<<it1->first<<"\t"<<it1->second<<endl;
////
////        }
//    }

    //copy_row(double *Y_gt, double *Y, int i, int num_rows, int num_cols)

}


#endif
