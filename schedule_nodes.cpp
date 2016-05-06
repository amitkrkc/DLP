void schedule_nodes(vector<vertex_vector>& J_all, const graph& g, const int& NUM_PROCESSORS)
{
	// clear J_all
	J_all.clear();

	vertex_vector J;

	vertex_set L;

	//initialize L
	for(tie(vi,vi_end)=vertices(g);vi!=vi_end;vi++)
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
			// if A is empty, break the for loop
                	if(!A.empty())
                	{
				// take te first vertex p and save it to J
                    		vertex_descriptor p=*(A.begin());
                    		J.push_back(p);

				// remove out-neighbors of p
                    		out_edge_iterator ei,ei_end;
				for(tie(ei,ei_end)=out_edges(p,g);ei!=ei_end;ei++)
					A.erase(target(*ei,g));            //average case complexity: linear in number of elements removed

				//remove in-neighbors of p
				in_edge_iterator ii,ii_end;
				for(tie(ii,ii_end)=in_edges(p,g);ii!=ii_end;ii++)
					A.erase(source(*ii,g));

				//remove p
				A.erase(p);
                	}
                	else
                    		break;
            	}

		// save J
		J_all.push_back(J);

		// remove the vertices in J from L
		for(int i=0;i<J.size();i++)
			L.erase(J[i]);
	}
}
