/*==========================================================
 * 
 * Author: ari pakman
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include <iostream>
#include <math.h>
#include <stdint.h>
#include "mex.h"

#include "HMC_BinarySampler.h"
#include "BinaryDistribution.h"
#include "MRF.h"

using namespace std;

//void _main();

/****************************/

static
void hmc_binary(

  double * samples, 
  double * loglikes, 
  double * M,
  double *r,
  int L,
  int P,
  int d, 
  double * last_Y,
  int * seed_ptr ){


  vector < vector<int> >  Ss;
  vector<double> log_likes; 

  vector <double>  Y(d);
  for (int i=0; i< d; i++)
    Y[i] = last_Y[i];

  
  MRF * mrf  = new MRF(d, M, r);


  HMC_BinarySampler sampler = HMC_BinarySampler(mrf);

  sampler.runSampler(L, P, Ss,log_likes, seed_ptr, Y);


  for (int i=0; i< d; i++)
    last_Y[i] = Y[i];


                         
  //copy the samples and loglikes to the variables to be returned
  for(int i=0; i < L; i++) {

    loglikes[i] = log_likes[i];
    
    for (int j=0; j < d; j++) {
      
      samples[d*i+j ] =  Ss[i][j];      
    }
  }




  return;
}


void mexFunction(
		 int          nlhs,
		 mxArray     *plhs[],
		 int          nrhs,
		 const mxArray *prhs[]
		 )
{
 
  
    // Input variables
    double *M;                /* 1xN input matrix */
    double *r;                /* 1xN input matrix */
    int L;
    int P;
    double *last_Y;
    int * seed_ptr = NULL;

    // Output variables
    double *samples;                   /* binary samples */
    double *loglikes;            /* log likelihood of each sample */


  /* Check for proper number of arguments */

  if (nrhs != 5 && nrhs != 6) {
    cout << nrhs << endl;
    mexErrMsgIdAndTxt("MATLAB:mexcpp:nargin", 
            "MEXCPP requires six or seven input arguments.");
  } else if (nlhs > 2) {
    mexErrMsgIdAndTxt("MATLAB:hmc_binary:nargout",
            "hmc_binary yields up to two output arguments.");
  }


  M = (double *) mxGetPr(prhs[0]);
  r = (double *) mxGetPr(prhs[1]);
  L = (int ) mxGetScalar(prhs[2]);
  P = (int ) mxGetScalar(prhs[3]);
  last_Y = (double *) mxGetPr(prhs[4]);
  if (nrhs == 6)
     seed_ptr = (int * ) mxGetPr(prhs[5]);




  /* check that number of cols in r is 1 */


  if(mxGetN(prhs[1])!=1) {
     mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input r must be a column vector.");
  }

  int d = mxGetM(prhs[1]);
  if(mxGetM(prhs[0])!=d || mxGetN(prhs[0])!=d) {
     mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","Input M must be a square matrix with same dim as r.");
  }


  // Create variables to return
  plhs[0] = mxCreateNumericMatrix(d, L, mxDOUBLE_CLASS, mxREAL);
  plhs[1] = mxCreateNumericMatrix(1, L, mxDOUBLE_CLASS, mxREAL);

  samples        = mxGetPr(plhs[0]);
  loglikes = mxGetPr(plhs[1]);

  hmc_binary(samples, loglikes, M,r,L,P, d, last_Y, seed_ptr);

  return;
}
