#ifndef SIMPLEX_H
#define SIMPLEX_H


#include "misc.h"

class simplex
{
    private:
        double c;
        vector_double v;
        double function_simplex(double nu) const;

    public:
        simplex(){}
        simplex(const vector_double& v1, const double& c1):v(v1), c(c1){}
        vector_double solve_quadratic_simplex(const options& opts);
};

// computes the function f:=\sum_i x_i-1 where x_i:=max(0, (v_i-nu)/c)
double simplex::function_simplex(double nu) const
{
	double f=-1.0;
	for(uint i=0;i<v.size();i++)
	{
		f+=std::max<double>(0.0, (v[i]-nu)/c);
	}
	return f;
}

// ******************************minimize the following problem************************************************************
//                  z=\arg\!\min_{z \in \Delta} \frac{1}{2}c_i z^\top z - v_i^\top z
// The strategy is to check the sign of c_i. If it is negative, the problem corresponds to
// minimizing a concave function over a probability simplex. If it is very small (==0), it reduces to a linear problem.
// The solution for (ci<0 or abs(ci)<1e-5) is obtained at one of the vertices of the simplex.
// When c_i>0, we solve it by using the KKT condition.

vector_double simplex::solve_quadratic_simplex(const options& opts)
{
    vector_double z(v.size(), 0.0);

    if(c<=0 || fabs(c)<1e-5)	//if c<=0, solution is obtained at one of the vertices
	{
		if(opts.verbose==2)
		{
            OUTPUT<<"\tc="<<c<<"\n\tminimize the concave funtion...solution is obtained at one of the vertices...\n";
		}

		uint minid=0;
		double minval=-v[0];
		for(uint i=1;i<v.size();i++)
		{
			double val=-v[i];
			if(val<minval)
			{
				minid=i;
				minval=val;
			}
		}

        z[minid]=1.0;
	}
	else
	{
		// minimize the convex program

		//compute max and min of v
		double max_v=v[0], min_v=v[0];
		for(uint i=1;i<v.size(); i++)
		{
			max_v=std::max<double>(max_v, v[i]);
			min_v=std::min<double>(min_v, v[i]);
		}

		//find the solution by bisection method
		double nu_1=max_v;
		double nu_0=min_v-c/v.size();

		double f1=function_simplex(nu_1);
		double f0=function_simplex(nu_0);
		if(opts.verbose==2)
			OUTPUT<<"nu_1:"<<nu_1<<"\tnu_0:"<<nu_0<<"\tf0:"<<f0<<"\tf1:"<<f1<<std::endl;

		double nu_m;

		while(1)
		{
			nu_m=(nu_0+nu_1)/2;
			double fm=function_simplex(nu_m);
			if(fabs(f1-f0)<opts.tol || fabs(fm)<opts.tol)
				break;
			else
			{
				if(f1*fm>0)
				{
					f1=fm;
					nu_1=nu_m;
				}
				else
				{
					f0=fm;
					nu_0=nu_m;
				}
			}

			if(opts.verbose==2)
				OUTPUT<<"\tnu_m:"<<nu_m<<"\tfm:"<<fm<<"\tf0:"<<f0<<"\tf1:"<<f1<<std::endl;
		}
		for(int i=0;i<z.size();i++)
		{
			z[i]=std::max<double>(0.0,(v[i]-nu_m)/c);
		}
	}
	return z;
}


#endif
