#ifndef REDCELL_H
#define REDCELL_H

#include <set>
#include <iostream>
#include <algorithm>
#include <vector>
#include <utility>      // for pair
#include <map>
#include "misc.h"

//using namespace std;          // it is NOT advised to use "using namespace" in header files

typedef std::pair<int,int> PAIR;                  //first int-->superpixel, second int-->area of the superpixel
typedef std::pair< std::vector<PAIR>, int> CELL;       // first--> pair of superpixel and area, second-->total area
typedef std::map<int, CELL> REDCELL;              // map is preferred to maintain the uniqueness of the red-cell label


// used for sorting the red cells or the superpixels in decreasing order of pixel area
template<class T>
struct greater_area_key
{
	bool operator() (const T& a, const T& b) const
	{
		return a.second>b.second;
	}
};


class Redcell
{
	private:
		void init(int *lut_data, const options& opts);		// private member function

	public:
        std::set<int> not_red_cells;     //stores the superpixels that do not belong to the red cells
		std::vector<CELL> redcells;		// vector that stores all red-cell superpixels

		Redcell(){}
		Redcell(int *lut_data, const options& opts)
		{
			redcells.clear();
            not_red_cells.clear();

			// calls init()
			init(lut_data, opts);
		}

		void display() const;

		int get_num_not_red_cells() const
		{
            return not_red_cells.size();
		}
};

void Redcell::display() const
{
    OUTPUT<<"***********REDCELLS*************"<<std::endl
        <<"num_not_red_cells:"<<get_num_not_red_cells()<<std::endl;
    OUTPUT<<"not_red_cells:";
    for(std::set<int>::const_iterator it=not_red_cells.begin(); it!=not_red_cells.end(); ++it)
        OUTPUT<<*it<<"\t";
    OUTPUT<<std::endl;

    int iter=0;
    for(std::vector<CELL>::const_iterator it=redcells.begin(); it!=redcells.end(); ++it)
    {
        ++iter;
        OUTPUT<<"Red cell:"<<iter<<"\t total area:"<<it->second<<std::endl;
        for(std::vector<PAIR>::const_iterator it1=it->first.begin(); it1!=it->first.end(); ++it1)
        {
            OUTPUT<<"\t("<<it1->first<<","<<it1->second<<")"<<std::endl;
        }
    }
}


// reads from the inputs lut_data and inserts into REDCELL. Afterwards, it copies the data to a vector<cell> and uses total_area of the cell to sort in decreasing order
// lut_data is a Nx3 integer matrix where each row corresponds to (sp, curr_label, pix_area)
// opts is a structure that contains all options/parameters
void Redcell::init(int *lut_data, const options& opts)
{
	if(opts.verbose>1)
        OUTPUT<<"Inside Redcell::init()"<<std::endl;

	REDCELL lut;

    // first we add all the nodes in the not-red-cells set
    // later, we will remove the nodes that belong to the red cells

    if(opts.verbose>1)
        OUTPUT<<"Adding "<< opts.num_nodes<<" nodes in not_red_cells set"<<std::endl;

    for(int i=0;i<opts.num_nodes;i++)
        not_red_cells.insert(i);



	for(int i=0;i<opts.lut_size; ++i)
	{
		// read the current entry from the lut_data
		int sp=lut_data[i]-1;                       //NOTE: 1 is subtracted for C convention
		int curr_label=lut_data[i+opts.lut_size];
		int pix_area=lut_data[i+2*opts.lut_size];

		if(opts.verbose>1)
            OUTPUT<<"sp:"<<sp<<"\t curr label:"<<curr_label<<"\t pix area:"<<pix_area<<std::endl;

		// mark the is_red flag of the sp superpixel as true
        not_red_cells.erase(sp);

		//insert the data to REDCELL lut
		REDCELL::iterator it=lut.find(curr_label);

        //if the curr_label does not exist in the map, create a new one
        if(it==lut.end())
        {
            if(opts.verbose>1)
                OUTPUT<<"\t does not exist in the REDCELL...creating a new one"<<std::endl;
            //create a new cell object temp
            std::vector<PAIR> temp;
            temp.push_back(std::make_pair(sp, pix_area));
            lut[curr_label]=std::make_pair(temp, pix_area);
        }
        //otherwise update the entry in the map
        else
        {
            if(opts.verbose>1)
                OUTPUT<<"\t exists in the REDCELL...updating"<<std::endl;
            it->second.second+=pix_area;
            it->second.first.push_back(std::make_pair(sp, pix_area));
        }
	}//end for i=0 to lut_size-1


	// since map does not sort the red cells based on their area,
    // we first convert it to a vector and then invoke the sort function

    if(opts.verbose>1)
        OUTPUT<<"Converting from REDCELL to vector<cell>"<<std::endl;

    for(REDCELL::iterator it=lut.begin(); it!=lut.end(); ++it)
    {
        redcells.push_back(it->second);
    }

    // sort the redcells based on the area
    std::sort(redcells.begin(), redcells.end(), greater_area_key<CELL>());


    // sort the superpixels of each red cell based on the area

    if(opts.verbose>1)
        OUTPUT<<"Sorting the superpixels based on their areas"<<std::endl;

    for(std::vector<CELL>::iterator it=redcells.begin(); it!=redcells.end(); ++it)
    {
        std::sort(it->first.begin(), it->first.end(), greater_area_key<PAIR>());
    }

    if(opts.verbose>1)
        OUTPUT<<"Returning from Redcell::init()"<<std::endl;
}

//void Redcell::compute_gt_label(double *Y_gt, double* Y, REDCELL lut, const options& opts)
//{
//
//    //pre-computing some factors
//    const double factor_red_cells=-0.25/num_red_cells;
//    const double factor_not_red_cells=1.0/num_red_cells/num_not_red_cells;
//
//	const int label_vector_size=1+opts.label_vector_size;
//
//    // copy Y_gt to zero initially
//    // later, each row of Y_gt will be modified depending on whether the super-pixel belongs to the red-cell or not
//    // we use memset because we will modify only few columns of the superpixels that belong to the red-cell
//    memset(Y_gt, 0.0, sizeof(double)*opts.num_nodes*label_vector_size);
//
//
//    vector<int> max_labels; // required later for not-red-cell superpixels
//
//    // initially make all labels available
//    set<int> available_labels;
//    for(int i=1;i<=label_vector_size;i++)
//        available_labels.insert(i);
//
//    // for each red cell (which is sorted based on the area)
//    for(vector<cell>::iterator it=redcells.begin(); it!=redcells.end(); it++)
//    {
//        if(VERBOSE>1)
//        {
//                OUTPUT<<it->total_area<<"\t superpixels:{";
//                for(vector< PAIR >::iterator it1=it->sp_area.begin(); it1!=it->sp_area.end(); it1++)
//                    OUTPUT<<" ("<<it1->first<<","<<it1->second<<") ";
//                OUTPUT<<"}"<<endl;
//
//        }
//
//        //find the maxid that has the largest area
//        vector< PAIR > sp_maxid;
//        for(vector< PAIR >::iterator it1=it->sp_area.begin(); it1!=it->sp_area.end(); it1++)
//        {
//            int sp=it1->first;
//            int area=it1->second;
//            int id=argmax(Y, sp, opts.num_nodes, label_vector_size);
//            sp_maxid.push_back(make_pair(area, id));
//        }
//
//        // sort the sp_maxid in decreasing order of area
//        std::sort(sp_maxid.begin(), sp_maxid.end(),greater_area_maxid());
//
//        int maxid=-1;
//        for(vector< PAIR >::iterator it1=sp_maxid.begin(); it1!=it->sp_maxid.end(); it1++)
//        {
//            maxid=it->second;
//
//            //check if maxid is available
//            set<int>::iterator it_av=available_labels.find(maxid);
//
//            // if the label is available
//            if(it_av!=available_labels.end())
//            {
//                break;
//            }
//        }
//
//
//
//        if(VERBOSE>1)
//            OUTPUT<<"\tmaxid:"<<maxid<<"\t maxarea:"<<maxarea<<endl;
//
//        // if maxid is not -1
//        if(maxid!=-1)           //sanity check
//        {
//            // set all the red pixels to a vector with -1/4R at the maxid-th index and zero everywhere where R is the number of superpixels in the red cells
//            for(set<int>::iterator it1=it->superpixels.begin(); it1!=it->superpixels.end(); it1++)
//            {
//                if(VERBOSE>1)
//                    OUTPUT<<"\t setting "<<*it1<<" row to e-vector at "<<maxid<<endl;
//                //set_row(Y_gt, *it1, maxid, num_nodes, label_vector_size);
//                //DISABLED: set_row because Y_gt has already been memset to 0 before the beginning of the for loop
//                // set_row also sets the row to zero. save extra loop,
//                int sp=*it1;
//                *(Y_gt+maxid*num_nodes+sp)=factor_red_cells;
//                max_labels.push_back(maxid);
//            }
//        }
//    }
//
//    // for the superpixels that are not in the red-cells
//    for(set<int>::iterator it=set_not_red_cells.begin(); it!=set_not_red_cells.end(); it++)
//    {
//        int not_sp=*it;
//        for(uint j=0;j<max_labels.size();j++)
//        {
//            int lbl=max_labels[j];
//            Y_gt[not_sp+num_nodes*lbl]=factor_not_red_cells;
//        }
//    }
//}
//
//
//
//
//class Redcell
//{
//	private:
//		vector<redcell> vv;
//		set<int> not_red_cells;
//	public:
//		Redcell(){}
//
//};


#endif
