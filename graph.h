#ifndef GRAPH_H
#define GRAPH_H

#include "simplex.h"
#include <boost/graph/graph_traits.hpp>     // for boost graph
#include <boost/graph/adjacency_list.hpp>   // for boost graph
#include "boost/tuple/tuple.hpp"            // for tie
#include <sys/time.h>                       // for timing different portions of the code
#include <set>                              // for set
#include <utility>                          // for pair

/////////////////////////////////////////////////////////////////////////////
struct node
{
	vector_double label;        //label to be assigned to each node
    double score;               // alpha*p_i-beta*(1-p_i). This is related to the score of a node to be a false positive

    // constructors
	node(){}
	node(const vector_double& l, const double& s)
	{
		label=l;
		score=s;
	}
	void display();
};

void node::display()
{
    OUTPUT<<"\tscore:"<<score<<std::endl;
    OUTPUT<<"\tlabel:";
    for(vector_double::iterator it=label.begin(); it!=label.end(); it++)
        OUTPUT<<*it<<"\t";
    OUTPUT<<std::endl;
}

//---------------------------------------------------
struct arc
{
	double weight;          // weight of an edge. It can be negative too
	arc(double w=0.0)
	{
		weight=w;
	}
	void display()
	{
        OUTPUT<<"\t weight:"<<weight<<std::endl;
	}
};


//------------------typedef-------------------------
typedef boost::adjacency_list<boost::listS,boost::vecS,boost::undirectedS,node,arc > graph_t;
typedef boost::graph_traits<graph_t>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<graph_t>::vertex_iterator vertex_iterator;
typedef boost::graph_traits<graph_t>::edge_iterator edge_iterator;
typedef boost::graph_traits<graph_t>::in_edge_iterator in_edge_iterator;
typedef boost::graph_traits<graph_t>::out_edge_iterator out_edge_iterator;
typedef boost::graph_traits<graph_t>::edge_descriptor edge_descriptor;
typedef std::vector<vertex_descriptor> vertex_vector;
//typedef std::list<vertex_descriptor> vertex_list;
typedef std::pair<int,int> edge_t;
typedef std::set<int> set_int;
typedef std::set<vertex_descriptor> vertex_set;

//---------------------------------------------------
class Graph
{
    private:
        graph_t* g;
        //vertex_vector vlist;     //pre-allocating vlist (for speed?)

    public:
        //constructor
        Graph()
        {
            g=new graph_t;
            //vlist.clear();
        }

        Graph(double *edgelist, int *nuclei, double *score, const options& opts)
        {
            g=new graph_t;
            add_node_and_edges(edgelist, nuclei, score, opts);
        }

        // destructor
        ~Graph()
        {
            delete g;
        }

        // adds nodes and edges to the graph
        // inputs: edge list (mx3 matrix), nuclei list (kx1 matrix), score list (nx1 matrix) and options,  where n=# nodes, m=#edges, k=#nuclei
        void add_node_and_edges(double *edgelist, int *nuclei, double *score, const options& opts);

        // assembles all label distributions to a label assignment matrix Y of size nxk where n=#nodes, k=label vector size
        void get_labels(double* Y);


        //schedules the nodes of the graph for parallel processing
        // input: number of processors (const)
        //        J_all-->the vector of vector of vertices. Inner vector contains the vertices that can be run in parallel.
        void schedule_nodes(std::vector<vertex_vector>& J_all, const int& NUM_PROCESSORS);


        // displays graph
        void display();

        // gathers labels at a node
        void gather_labels_node(vector_double& v, double& c, vertex_descriptor p,  const options& opts);

        // solves the non-convex quadratic problem at the p-th node
        // returns the resulting label distribution as vector<double>
        // during training phase, y_gt will contain inforamtion from the red cells
        // during testing phase, y_gt will be a vetor of zeros
        vector_double solve_dc_nodewise(vertex_descriptor p,const vector_double& y_gt, const options& opts);


        // computes energy of the graph
        double compute_energy();


        // set the label distribution of the node p
        void set_label(vertex_descriptor p, const vector_double& z);

        //computes argmax of label distribution of the node p
        int argmax(const vertex_descriptor& p)
        {
            vector_double::iterator it=std::max_element((*g)[p].label.begin(), (*g)[p].label.end());
            return static_cast<int>(it-(*g)[p].label.begin());
        }

        int argmax(const int& p)
        {
            vector_double::iterator it=std::max_element((*g)[p].label.begin(), (*g)[p].label.end());
            return static_cast<int>(it-(*g)[p].label.begin());
        }

        double compute_label_difference(int p, const vector_double& z)
        {
            return norm2((*g)[p].label, z);
        }

        //void perform_label_propagation(const options& opts);
};
///////////////////////////////////////////////////////////////
// performs label propagation
//void Graph::perform_label_propagation(const options& opts)
//{
//
//}

//////////////////////////////////////////////////////////////

//sets the label distribution of the p-th node
void Graph::set_label(vertex_descriptor p, const vector_double& z)
{
    assert(z.size()==(*g)[p].label.size());
    (*g)[p].label=z;
}

//////////////////////////////////////////////////////////////

double Graph::compute_energy()
{
    // total energy f, initialized to zero
	double f=0.0;

    // label vector size
    const uint label_vector_size=(*g)[0].label.size();

    // false positive label
	vector_double l_0(label_vector_size, 0.0);
	l_0[0]=1.0;

    // iterate through all vertices
    vertex_iterator vi, vi_end;
	for(boost::tie(vi, vi_end)=vertices(*g); vi!=vi_end; vi++)
	{
        vertex_descriptor p=*vi;

        // label of the p-th vertex
        vector_double y_p=(*g)[p].label;

        // L2 norm of y_p-L_0
        double term2=0.5*(*g)[p].score*norm2(y_p, l_0);

        // label dissimilarity between the edges
        double term3=0.0;
        out_edge_iterator ei, ei_end;
        for(boost::tie(ei, ei_end)=out_edges(p,*g); ei!=ei_end; ei++)
        {
            vertex_descriptor j=target(*ei,*g);     //node j
            vector_double y_j=(*g)[j].label;        // label of the j-th node
            double w_pj=(*g)[*ei].weight;           // weight of the (p,j) edge
            term3+=w_pj*norm2(y_p, y_j);            // label discrepancy between y_p and y_j
        }

        term3*=0.5;

        f+=term2+term3;
	}
	return f;
}

// gathers the labels of all the neighboring nodes of the given node p
void Graph::gather_labels_node(vector_double& v_p, double& c_p, vertex_descriptor p, const options& opts)
{

    assert(v_p.size()==(*g)[p].label.size());

    // initialize v_p
    std::fill(v_p.begin(), v_p.end(), 0.0);
    v_p[0]+=(*g)[p].score;                 // add the score

	// gather labels from the neighbors
	// initialize c_p
	c_p=(*g)[p].score;

	//edges of p
	out_edge_iterator ei,ei_end;
	for(boost::tie(ei,ei_end)=out_edges(p,*g);ei!=ei_end;ei++)
	{
		vertex_descriptor t=target(*ei,*g);
		vector_double y_j=(*g)[t].label;
		double w_pj=(*g)[*ei].weight;

		if (opts.verbose>1)     //2--> full verbose
		{
			OUTPUT<<"source:"<<p<<"\t target:"<<t<<"\t weight:"<<w_pj<<std::endl;
			for(vector_double::iterator it=y_j.begin(); it!=y_j.end(); it++)
                OUTPUT<<*it<<"\t";
            OUTPUT<<std::endl;
		}

		c_p+=w_pj;      // add w_pj

		// gather labels
		for(uint i=0;i<y_j.size();i++)
            v_p[i]+=w_pj*y_j[i];
	}
}

/////////////////////////////////////////////////////////////

// displays the graph
void Graph::display()
{
    OUTPUT<<"******************************\n"
            <<"Num. nodes:"<<num_vertices(*g)<<std::endl
            <<"Num. edges:"<<num_edges(*g)<<std::endl
            <<"label vector size:"<<(*g)[0].label.size()<<std::endl
            <<"--------------------------------"<<std::endl;

    // display labels of all the nodes
    vertex_iterator vi, vi_end;
    for(boost::tie(vi, vi_end)=vertices(*g); vi!=vi_end; vi++)
    {
        OUTPUT<<"Node:"<<*vi<<std::endl;
        (*g)[*vi].display();
    }

    OUTPUT<<"-------------------------"<<std::endl;

    // display edges
    out_edge_iterator ei, ei_end;
    for(boost::tie(vi, vi_end)=vertices(*g); vi!=vi_end; vi++)
    {
        for(boost::tie(ei, ei_end)=out_edges(*vi, *g); ei!=ei_end; ei++)
        {
            vertex_descriptor s=source(*ei, *g), t=target(*ei, *g);
            OUTPUT<<"Edge: ("<<s<<","<<t<<")";
            (*g)[*ei].display();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

// solve label propagation at node p
vector_double Graph::solve_dc_nodewise(vertex_descriptor p,const vector_double& y_gt, const options& opts)
{
    assert(y_gt.size()==opts.label_vector_size+1);      //+1 because of the false positive label

    if(opts.verbose>1)
        OUTPUT<<"Inside solve_dc_nodewise function\n";

    // initialize v_p
    vector_double v_p(y_gt.size(),0.0);
    double c_p=0.0;

    // gather labels from the neighbors
    // remember that this does not incorporate the red cells yet
    gather_labels_node(v_p, c_p, p, opts);

    // during training phase y_gt will be obtained from the red-cell
    // during testing phase, y_gt will be a vector of zeros
    for(int i=0;i<y_gt.size();i++)
    {
        v_p[i]+=y_gt[i];
    }

	if(opts.verbose>1)
	{
		OUTPUT<<"-----------------------------------------\n";
		OUTPUT<<"c_p:"<<c_p<<std::endl;
		for(vector_double::iterator it=v_p.begin(); it!=v_p.end(); it++)
            OUTPUT<<*it<<"\t";
//		display<double>(&v[0], v.size());
		OUTPUT<<"\n-----------------------------------------\n";
	}

	// solve the quadratic problem at node p
	simplex s(v_p, c_p);                                // create a simplex object
	vector_double z=s.solve_quadratic_simplex(opts);	// z contains the updated label of node p

	if(opts.verbose>1)
	{
		OUTPUT<<"Returning from optimization\t optimal label:";
		for(vector_double::iterator it=z.begin(); it!=z.end(); it++)
            OUTPUT<<*it<<"\t";
        OUTPUT<<std::endl;
//		display<double>(&z[0], z.size());
	}
	return z;
}

/////////////////////////////////////////////////////////////////////////
void Graph::get_labels(double* Y)
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

///////////////////////////////////////////////////////////////////////

void Graph::schedule_nodes(std::vector<vertex_vector>& J_all,const int& NUM_PROCESSORS)
{
	vertex_vector J;

	vertex_set L;

	//initialize L
	vertex_iterator vi, vi_end;
	for(boost::tie(vi,vi_end)=vertices(*g);vi!=vi_end;vi++)
	{
		L.insert(*vi);  //average case complexity of insertion: constant, worst case comlexity: linear in container size
	}

	while(!L.empty())
    {
        // clear J
		J.clear();

		//initalize with L
		vertex_set A(L);

		// try to find NUM_PROCESSORS nodes for parallel processing from A
		for(int i=0;i<NUM_PROCESSORS;i++)
        {
			// if A is not empty, pick a node
            if(!A.empty())
            {
				// take te first vertex p and save it to J
                vertex_descriptor p=*(A.begin());
                J.push_back(p);

				// remove out-neighbors of p
                out_edge_iterator ei,ei_end;
				for(tie(ei,ei_end)=out_edges(p,*g);ei!=ei_end;ei++)
				{
					A.erase(target(*ei,*g));            //average case complexity: linear in number of elements removed
                }

				//remove in-neighbors of p
				in_edge_iterator ii,ii_end;
				for(tie(ii,ii_end)=in_edges(p,*g);ii!=ii_end;ii++)
                {
					A.erase(source(*ii,*g));
                }
				//remove p
				A.erase(p);
            }
            else // if A is empty, break the for loop
            {
                break;
            }
        }

		// save J
		J_all.push_back(J);

		// remove the vertices in J from L
		for(int i=0;i<J.size();i++)
		{
			L.erase(J[i]);
        }
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

// add node(s) and edges from the input lists
void Graph::add_node_and_edges(double *edgelist, int *nuclei, double *score, const options& opts)
{
	if(opts.verbose>1)
		OUTPUT<<"\n\tCreating graph...\n";


    vertex_vector vlist(opts.num_nodes);


    // add all nuclei to the set
    // nuclei_list is a num_nodes x 1 array of booleans to store if a super-pixel is a nucleus or not.
    // it offers a constant time look-up later

    if(opts.verbose>1)
        OUTPUT<<"\tAdding nuclei...";
    std::vector<bool> nuclei_list(opts.num_nodes, false);
    for(int k=0;k<opts.num_nuclei;k++)
    {
        int q=nuclei[k]-1;      //note the -1 because of C convention
        nuclei_list[q]=true;
    }

    if(opts.verbose>1)
        OUTPUT<<"done"<<std::endl<<"\n\tAdding vertices...";

    	// add nodes to the graph
	for(int k=0;k<opts.num_nodes;k++)
	{
        	// initialize the label vector
		vector_double label(opts.label_vector_size+1,0.0);          //+1 because of the background label
		if(nuclei_list[k+1]==true)//is_nuclei(nuclei_list, k))       //should be efficient than the previous version. Compared to earlier set based solution, this solution is O(1). is_nuclei() is defined in misc.h
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
    if(opts.verbose>1)
        OUTPUT<<"done...added "<<num_vertices(*g)<<" vertices"<<std::endl<<"\t Adding edges...";

    //add edges
    double *srclist=edgelist;
    double *destlist=edgelist+opts.num_edges;
    double *weightlist=edgelist+2*opts.num_edges;

	for(int k=0;k<opts.num_edges;k++)
	{
		int s=(int)(*srclist++)-1;		//-1 because of c convention
		int t=(int)(*destlist++)-1;
		double w=*weightlist++;
		add_edge(vlist[s],vlist[t],arc(w),*g);
	}
	if(opts.verbose>1)
        OUTPUT<<"done"<<std::endl;
}



#endif
