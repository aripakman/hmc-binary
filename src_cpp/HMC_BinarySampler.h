/* 
 * File:   HMC_BinarySampler.h
 * Author: aripakman
 *
 * Created on May 8, 2013, 10:37 PM
 */

#ifndef HMC_BINARYSAMPLER_H
#define	HMC_BINARYSAMPLER_H

#include <vector>
#include "BinaryDistribution.h"


using namespace std;


class HMC_BinarySampler {
    
public:
   
    HMC_BinarySampler(const BinaryDistribution * b ): bin_dist(b) {}     // b is passed as a pointer because an object of type BinaryDistribution 
                                                                         // cannot be created, since it contains pure virtual functions 
    
    void runSampler(const int & L, const int & P, vector<vector<int>>  &Ss, vector<double> & log_likes, int * seed_ptr,  vector<double> & Y);

    int get_dim() {return bin_dist->getDimension();}
    
private:
    
    double _getHitTime(const double & x, const double & v );
    
    const BinaryDistribution * bin_dist;

};

#endif	/* HMC_BINARYSAMPLER_H */

