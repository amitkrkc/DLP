#ifndef DLP_H
#define DLP_H

#include "graph.h"
#include "redcell.h"

class DLP:public Redcell, public Graph
{
    private:
//        Redcell redcells;
//        Graph g;
        void compute_gt_label(double *Y_gt, const options& opts);
    public:
        DLP(){}
        DLP(double *edgelist, int *nuclei, double *score, int *lut_data, const options& opts):Graph(edgelist, nuclei, score, opts), Redcell(lut_data, opts){}
        void perform_label_propagation(vector_double& energy_evo, const options& opts);
        //void get_labels(double *Y);
};
///////////////////////////////////////////////////////

//void DLP::get_labels(double *Y)
//{
//    g.get_labels(Y);
//}

///////////////////////////////////////////////////////

void DLP::perform_label_propagation(vector_double& fevo, const options& opts)
{
    if(opts.verbose)
    {
        OUTPUT<<"Inside DLP::perform_label_propagation function"<<std::endl;
        opts.display();
    }

    struct timeval t_begin, t_end;      // for measuring time taken by different steps

    // schedule nodes
    std::vector<vertex_vector> J_all;
    const int num_proc=opts.use_parallel_processing?opts.num_processors:1;
    if(opts.verbose>1)
    {
        OUTPUT<<"Scheduling the "<<opts.num_nodes<<" nodes on "<<num_proc<<" parallel processors"<<std::endl;
    }

    gettimeofday(&t_begin,NULL);

    schedule_nodes(J_all, num_proc);

    gettimeofday(&t_end,NULL);

    if(opts.verbose>1)
    {
        double elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
        OUTPUT<<"Scheduling took "<< elapsed_ms<<" milliseconds"<<std::endl;
    }

    double *Y_gt=new double[opts.num_nodes*(opts.label_vector_size+1)];


    vertex_iterator vi, vi_end;
    double max_change=-100000000000; //initialized to a very small number
    double fun_change=100000000;        //initialized to a very big number

    // main loop
    int iter=0;

    fevo.clear();
    while(1)
    {
        iter++;
        gettimeofday(&t_begin,NULL);

        compute_gt_label(Y_gt, opts);       // compute the gt label using red cell and current labels

        gettimeofday(&t_end,NULL);
        double elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;
        if(opts.verbose>1)
            OUTPUT<<"\t compute_gt_labels took "<<elapsed_ms<< " milliseconds"<<std::endl;

        gettimeofday(&t_begin,NULL);
        for(std::vector<vertex_vector>::iterator it=J_all.begin(); it!=J_all.end(); ++it)
        {
            for(vertex_vector::iterator it1=it->begin(); it1!=it->end(); ++it1)
            {
                vertex_descriptor p=*it1;

                vector_double y_gt(opts.label_vector_size, 0.0);    //default is zero
                if(opts.is_train_phase)
                    get_row(y_gt, static_cast<int>(p), Y_gt, opts);

                vector_double z_p=solve_dc_nodewise(p, y_gt, opts);

                max_change=std::max<double>(max_change, compute_label_difference(p,z_p));


                g.set_label(p, z_p);
            }

        }

        gettimeofday(&t_end,NULL);
        elapsed_ms=(t_end.tv_sec-t_begin.tv_sec)*1000.0+(t_end.tv_usec-t_begin.tv_usec)/1000.0;

        double eval=compute_energy();

		if(fevo.size()>1)
	        	fun_change=fabs(eval-*fevo.rbegin());
		fevo.push_back(eval);

        if(opts.verbose>0)
            OUTPUT<<"iteration:"<<iter<<"\tgraph energy:"<<eval<<"\t max change:"<<max_change<<"\t fun change:"<<fun_change<<"\ttime taken:"<<elapsed_ms<<std::endl;

		if(iter>opts.maxIter || fun_change<opts.tol || max_change<opts.tol)
			break;
    }


	// free the memory
	delete[] Y_gt;

    if(opts.verbose>0)
        OUTPUT<<"returning from perform_label_propagation"<<std::endl;

	return fevo;
}
////////////////////////////////////////////

void DLP::compute_gt_label(double *Y_gt, const options& opts)
{

    //pre-computing some factors
    const double factor_red_cells=-0.25/redcells.num_red_cells;
    const double factor_not_red_cells=1.0/redcells.num_red_cells/redcells.num_not_red_cells;

	const int label_vector_size=1+opts.label_vector_size;   //+1 because of the false positive label

    // copy Y_gt to zero initially
    // later, each row of Y_gt will be modified depending on whether the super-pixel belongs to the red-cell or not
    // we use memset because we will modify only few columns of the superpixels that belong to the red-cell
    memset(Y_gt, 0.0, sizeof(double)*opts.num_nodes*label_vector_size);


    vector<int> max_labels; // required later for not-red-cell superpixels

    // initially make all labels available
    // note that i is initialized to 1 to exclude the false positive label (0)
    std::set<int> available_labels;
    for(int i=1;i<=label_vector_size;++i)
        available_labels.insert(i);

    // for each red cell (which is sorted based on the area)
    for(vector<CELL>::iterator it=redcells.begin(); it!=redcells.end(); ++it)
    {

        // traverse the (sorted) superpixel-area vector
        int maxid=-1;
        for(vector< PAIR >::iterator it1=it->first.begin(); it1!=it->first.end(); ++it1)
        {
            // compute the argmax of the label of the superpixel it1->first
            maxid=g.argmax(it1->first);

            if(opts.verbose>1)
                OUTPUT<<"\t superpixel:"<<it1->first<<"\t area:"<<it1->second<<"\t maxid:"<<maxid<<std::endl;

            // check if maxid is available
            std::set<int>::iterator it_av=available_labels.find(maxid);
            if(it_av!=available_labels.end())
            {
                if(opts.verbose>1)
                    OUTPUT<<"\t label:"<<maxid<<" is available...breaking..."<<std::endl;

                break;
            }
        }

        // if maxid is not -1
        if(maxid!=-1)           //sanity check
        {
            max_labels.push_back(maxid);
            // set all the red pixels to a vector with -1/4R at the maxid-th index and zero
            // everywhere where R is the number of superpixels in the red cells
            for(vector< PAIR >::iterator it1=it->first.begin(); it1!=it->first.end(); ++it1)
            {
                int sp=it1->first;

                if(opts.verbose>1)
                    OUTPUT<<"\t setting "<<sp<<" row to e-vector at "<<maxid<<std::endl;

                *(Y_gt+maxid*opts.num_nodes+sp)=factor_red_cells;       // we put factor_red_cells instead of 1 because of the way the loss function is defined
                                                                        // loss function is defined as \frac{1}{|R|} \sum_{i in R} y_gt_i^\top y_i
            }

        }
    }

    // for the superpixels that are not in the red-cells
    for(set<int>::iterator it=set_not_red_cells.begin(); it!=set_not_red_cells.end(); ++it)
    {
        int not_sp=*it;
        for(uint j=0;j<max_labels.size();j++)
        {
            int lbl=max_labels[j];
            Y_gt[not_sp+num_nodes*lbl]=factor_not_red_cells;
        }
    }
}



#endif
