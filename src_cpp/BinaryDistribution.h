/* 
 * File:   BinaryDistribution.h
 * Author: aripakman
 *
 * Created on May 8, 2013, 9:57 PM
 */

#include<vector>

#ifndef BINARYDISTRIBUTION_H
#define	BINARYDISTRIBUTION_H


using namespace std;

class BinaryDistribution {
public:
    BinaryDistribution();
    BinaryDistribution(const BinaryDistribution& orig);
    virtual ~BinaryDistribution();
    
    virtual double getLogLike( vector<int> & ) const = 0;
    virtual double getLogLikeDifference(vector<int> &, int &) const = 0;
    virtual int getDimension() const {return d;}
    
protected:
    int d;    // dimension of the binary vector
};

#endif	/* BINARYDISTRIBUTION_H */

