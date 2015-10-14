/* 
 * File:   Ising2D.h
 * Author: aripakman
 *
 * Created on October 12, 2015
 */

#ifndef MRF_H
#define	MRF_H

#include <vector>
#include <array>
#include "BinaryDistribution.h"


using namespace std;



class MRF: public BinaryDistribution {
public:

	MRF(int d2, double *M2, double *r2) : d(d2), M(M2), r(r2) 	{}
	//MRF(double beta, int d, double *M, double *r);
    double getLogLike( vector<int> & ) const;
    double getLogLikeDifference(vector<int> &, int &) const;
    int getDimension() const {return d;}
    
    
private:    
    int d;
    double *M;
    double *r;    
};

#endif

