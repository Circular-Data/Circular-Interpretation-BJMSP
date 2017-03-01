/*
----------------------------------------------------------
VenterModeCircular.cpp
Calculate the Highest Posterior density for circular variables (HPD, as in Venter (1967), to either
estimate the mode or obtain the shortest credible interval, for example.

----------------------------------------------------------
*/

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// [[Rcpp::depends(BH)]]
#include <boost/math/special_functions/bessel.hpp>

#include <iostream>
#include <math.h>
#include <typeinfo>


using namespace Rcpp;
using namespace std;
using namespace arma;




// [[Rcpp::export]]
double hmodeC(NumericVector x, double cip) {
  /* FUNCTION hmode -------------------------------------------
  Estimate the mode by finding the highest posterior density interval.

  x:      Sample from which to estimate the mode.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.

  Returns: An scalar containing the estimate of the mode.
  ------------------------------------------------------------ */

  int n, cil, chiv;
  double ln, M;
  const double pi = boost::math::constants::pi<double>();
  
  n = x.size();
  NumericVector sx = clone(x);
  NumericVector sx2 = clone(x)+(2*pi);
  std::vector<double> SX;
  SX.reserve( x.size() + x.size() ); // preallocate memory
  SX.insert( SX.end(), sx.begin(), sx.end() );
  SX.insert( SX.end(), sx2.begin(), sx2.end() );
  std::sort(SX.begin(), SX.end()); 
  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Size of the currently smallest interval.
  ln = SX[cil]-SX[0];

  for (int i=0; i < (n); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (SX[i+cil]-SX[i])) {
      ln = (SX[i+cil]-SX[i]);
      chiv = i;
    }
  }

  M = (fmod(SX[chiv+cil],(2*pi))+SX[chiv])/2;

  return M;
}



// Find the highest density interval.

// [[Rcpp::export]]
NumericVector hmodeciC(NumericVector x, double cip) {
  /* FUNCTION hmodeci -----------------------------------------
  Find the highest posterior density interval.

  x:      Sample from which to estimate the interval.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.

  Returns: An vector of length 2 containing
           lower and upper bound of the interval.
  ------------------------------------------------------------ */

  int n, cil, chiv;
  double ln;
  const double pi = boost::math::constants::pi<double>();
  
  n = x.size();
  NumericVector sx = clone(x);
  NumericVector sx2 = clone(x)+(2*pi);
  std::vector<double> SX;
  SX.reserve( x.size() + x.size() ); // preallocate memory
  SX.insert( SX.end(), sx.begin(), sx.end() );
  SX.insert( SX.end(), sx2.begin(), sx2.end() );
  std::sort(SX.begin(), SX.end());
  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);

  // Will be the minimal value of the smallest interval.
  chiv = 0;

  // Length of the currently smallest interval.
  ln = SX[cil]-SX[0];

  for (int i=0; i < (n); i++) {

    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (SX[i+cil]-SX[i])) {
      ln = (SX[i+cil]-SX[i]);
      chiv = i;
    }
  }

  NumericVector M(2);
  M[0] = SX[chiv];
  M[1] = fmod(SX[chiv+cil],(2*pi));

  return M;
}
