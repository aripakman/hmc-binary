/* 
 * File:   HMC_BinarySampler.cpp
 * Author: aripakman
 * 
 * Created on May 8, 2013, 10:37 PM
 */

#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <random>
#include <iostream>
#include <vector>
#include <utility>  //for the pair templated class
#include "BinaryDistribution.h"
#include "HMC_BinarySampler.h"

using namespace std;

void HMC_BinarySampler::runSampler(const int & L, const int & P, vector<vector<int> >  &Ss, vector<double> & log_likes, int * seed_ptr, vector<double> & Y ){

    int d = bin_dist->getDimension();

    int seed;
    if (seed_ptr == NULL)
        seed = rand();
    else 
        seed = *seed_ptr;
    
    default_random_engine e(seed);
    normal_distribution<> normal(0,1);
   
    typedef pair<int,double> Hpair;    
 
    Ss.clear();
    log_likes.clear();


    vector<double> V(d);    // velocities
    vector<int> S(d);       //signs of Y    

    // sample initial positions
    for (int j=0; j < d; j++ ){                        
        S[j] = (Y[j] < 0) ? -1 : 1;
        
    }
   
    // store initial values
    // Ss.push_back(S);
    // log_likes.push_back(bin_dist->getLogLike(S));
    
    for (int i=0; i < L; i++ ){
        
    // sample initial velocities
        for (int j=0; j < d; j++ ){                        
            V[j] = normal(e);
        }
                
    // compute initial hit times       
        
        vector<Hpair> hit_times(d);
        
        for (int j=0; j < d; j++ ){                        
            hit_times[j] = make_pair(j,_getHitTime( Y[j], V[j] ));
        }
    
        // sort the coordinates according to their hit_times
        sort(hit_times.begin(), hit_times.end(), [] (const Hpair &pair1, const Hpair & pair2 ) {return pair1.second < pair2.second; }   );
       
        
    // move the particle 
      
    // the first T =pi time causes every coordinate to reach zero 
        for (int j= 0; j < d; j++){        
            
                int c = hit_times[j].first;
                double t = hit_times[j].second;
                V[c] = V[c]* cos(t) - Y[c]*sin(t);  
                
                double V2new = V[c]*V[c] - 2*S[c]*bin_dist->getLogLikeDifference(S ,c);                    
                if (V2new > 0){
                    V[c] = -S[c]*sqrt(V2new);
                    S[c] = -S[c];
                } else {
                    V[c] = - V[c];                    
                }
         }
        
        
     // The next P-1 cycles of T=pi have known hit times and hit velocities   
        
        for (int p = 1; p < P; p++){
            for (int j= 0; j < d; j++){
                
                int c = hit_times[j].first;                                
                
                double V2new = V[c]*V[c] - 2*S[c]*bin_dist->getLogLikeDifference(S ,c);
                if (V2new > 0){
                    V[c] = -S[c]*sqrt(V2new);
                    S[c] = -S[c];
                } 
                // if the particle does not cross there are here two sign inversions for V[c]. 
                // the first one is from moving t = pi, the second for being reflected.

                        
                    
            }  // for j                                  
        } // for p
        


     // At this point, all the coordinates have moved a time hit_times[j].second + (P-1)*pi   
     // and need to move 1.5*pi - hit_times[j].second. 
     
       
            for (int j= 0; j < d; j++){
                
                int c = hit_times[j].first;                                
                double mt; // final_move
                if (hit_times[j].second < M_PI/2){
                    
                    // move t=pi    
                    double V2new = V[c]*V[c] - 2*S[c]*bin_dist->getLogLikeDifference(S ,c);
                    if (V2new > 0){
                        V[c] = -S[c]*sqrt(V2new);
                        S[c] = -S[c];
                    } 

                    mt = M_PI/2 - hit_times[j].second;
                }
                else {
                   mt= 3*M_PI/2 - hit_times[j].second;}

                Y[c] = V[c]*sin(mt);
                    
            }  // for j = 0 to d-1                                  

        
        // store the signs of X and the log_likelihood of the state
       
        Ss.push_back(S);
        log_likes.push_back(bin_dist->getLogLike(S));
        
    } // for i =1 to L -1

}


double HMC_BinarySampler::_getHitTime(const double & y, const double & v ){
    
    double phi = atan2(y, v);           // - pi < phi < +pi
    double t1;   // = -phi + M_PI/2;           // first time to hit x=0
                                         // - pi/2 < phi < +3pi/2
    if ( phi> 0) {t1 = M_PI -phi ;}
    else {t1 = -phi;}
        
    
    return t1;    
}

