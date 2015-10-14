/* 
 * File:   Ising2D.cpp
 * Author: aripakman
 * 
 * Created on May 8, 2013, 10:18 PM
 */


#include <vector>
#include <iostream>
#include "MRF.h"


using namespace std;




/*MRF::MRF(double beta, int d, double *M, double *r): beta(beta), d(d), M(M), r(r) {} */


double MRF::getLogLike( vector<int> & S)  const {
    double ll=0;
    for (int i = 0; i<d; i++){
        ll  -= r[i]*S[i];

        for (int j=0; j<i; j++){
            ll -= M[i*d + j]*S[j]*S[i];
        }
    }                   
    return ll;
    }
    

double MRF::getLogLikeDifference(vector<int> & S, int &i) const {

    double diff = -2*r[i];

    for (int j=0; j< d; j++){
        if (j==i) continue;
        diff -= 2*M[i*d + j]*S[j];
    }
    return diff;
    

}
