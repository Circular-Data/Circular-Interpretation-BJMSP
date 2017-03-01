// [[Rcpp::depends(Rcpp, RcppArmadillo, BH)]]

#include <RcppArmadillo.h>
#include <boost/math/special_functions/bessel.hpp>
#include <math.h>
#include <iostream>
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double hmodeC(arma::vec x, double cip) {
  /* FUNCTION hmode -------------------------------------------
  Estimate the mode by finding the highest posterior density interval
  for a circular variable.
  
  x:      Sample from which to estimate the mode.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.
  
  Returns: An scalar containing the estimate of the mode.
  ------------------------------------------------------------ */
  
  int n, cil, chiv;
  double ln, M;
  const double pi = boost::math::constants::pi<double>();
  
  n = x.size();
  arma::mat sx = x;
  arma::mat sx2 = x+(2*pi);
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
arma::vec hmodeciC(arma::vec x, double cip) {
  /* FUNCTION hmodeci -----------------------------------------
  Find the highest posterior density interval for a circular variable.
  
  x:      Sample from which to estimate the interval.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.
  
  Returns: An vector of length 2 containing
  lower and upper bound of the interval.
  ------------------------------------------------------------ */
  
  int n, cil, chiv;
  double ln;
  const double pi = boost::math::constants::pi<double>();
  
  n = x.size();
  arma::vec sx = x;
  arma::vec sx2 = x+(2*pi);
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
  
  arma::vec M(2);
  M[0] = SX[chiv];
  M[1] = fmod(SX[chiv+cil],(2*pi));
  
  return M;
}

// [[Rcpp::export]]
double hmode(arma::vec x, double cip) {
  /* FUNCTION hmode -------------------------------------------
  Estimate the mode by finding the highest posterior density interval.
  
  x:      Sample from which to estimate the mode.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.
  
  Returns: An scalar containing the estimate of the mode.
  ------------------------------------------------------------ */
  
  int n, cil, chiv;
  double ln, M;
  
  n = x.size();
  arma::vec sx = x;
  std::sort(sx.begin(), sx.end());
  
  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);
  
  // Will be the minimal value of the smallest interval.
  chiv = 0;
  
  // Size of the currently smallest interval.
  ln = sx[cil]-sx[0];
  
  for (int i=0; i < (n-cil); i++) {
    
    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }
  
  M = (sx[chiv+cil]+sx[chiv])/2;
  
  return M;
}

// Find the highest density interval.

// [[Rcpp::export]]
arma::vec hmodeci(arma::mat x, double cip) {
  /* FUNCTION hmodeci -----------------------------------------
  Find the highest posterior density interval.
  
  x:      Sample from which to estimate the interval.
  cip:    Bandwith for the algorithm, ranging from 0 to 1.
  
  Returns: An vector of length 2 containing
  lower and upper bound of the interval.
  ------------------------------------------------------------ */
  
  int n, cil, chiv;
  double ln;
  
  n = x.size();
  arma::mat sx = x;
  std::sort(sx.begin(), sx.end());
  
  // The number of values within the
  // (cip*100)% Confidence Interval
  cil = trunc(cip*n);
  
  // Will be the minimal value of the smallest interval.
  chiv = 0;
  
  // Length of the currently smallest interval.
  ln = sx[cil]-sx[0];
  
  for (int i=0; i < (n-cil); i++) {
    
    // If the smallest interval so far is larger than the
    // current, set the current as the new smallest interval.
    if (ln > (sx[i+cil]-sx[i])) {
      ln = (sx[i+cil]-sx[i]);
      chiv = i;
    }
  }
  
  arma::vec M(2);
  M[0] = sx[chiv];
  M[1] = sx[chiv+cil];
  
  return M;
}


// [[Rcpp::export]]
//Eigenvalues for a dense square matrix
arma::vec Eigenvalues(arma::mat X) {
  
  arma::cx_vec eigval;
  arma::eig_gen(eigval, X);
  
  return arma::conv_to<arma::vec>::from(eigval);
}

// [[Rcpp::export]]
//Eigenvectors for dense general square nonsymmetric matrix
arma::mat Eigenvectors(arma::mat X) {
  
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, X);
  
  arma::mat eigvc = arma::conv_to<arma::mat>::from(eigvec);
  return  -1 * eigvc;
}

// [[Rcpp::export]]
//Eigenvectors for a dense symmetric matrix
arma::mat EigenvectorsSym(arma::mat X) {
  
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, X, "dc");
  
  arma::mat eigvc = arma::conv_to<arma::mat>::from(eigvec);
  return  -1 * eigvc;
}

// [[Rcpp::export]]
//Check if a matrix is identity
//Note that this function only works for square matrices
int IsIdentity(arma::mat X){
  //count rows and columns in X
  int ncols = X.n_cols;
  int nrows = X.n_rows;
  
  //create a testvector 
  
  arma::vec testI(ncols*nrows);
  
  //set a count
  int count = 0;
  
  //start checking
  for(int iii = 0; iii < ncols; ++iii){
    for(int jjj = 0; jjj < nrows; ++jjj){
      int testsim = X(jjj,iii);
      int testsim2 = X(iii,jjj);
      if(((testsim == 1) & (iii == jjj)) | ((testsim2 == 0) & (iii != jjj))){
        testI(count) = 0;
      }else{
        testI(count) = 1;
      }
      count = count + 1;
    }
  }
  //return sum of the test vector 
  //if sum of this vector is zero, the matrix is an identity matrix
  return sum(testI);
}

// [[Rcpp::export]]
//Sample from a multivariate normal distribution
//Note that sigma has to be a square matrix
arma::mat mvrnormArmaEigen(int n, arma::vec mu, arma::mat sigma){
  int ncols = sigma.n_cols;
  arma::mat Y(n, ncols);
  for(int iii=0; iii < ncols; ++iii){
    Y.col(iii) = as<arma::vec>(Rcpp::rnorm(n));
  }
  
  int testsym = 0;
  
  //test if matrix is idenity matrix
  if(IsIdentity(sigma) == 0){
    testsym = 1;
  }
  
  //test is matrix is symmetrical
  arma::mat sigmat = sigma.t();
  for(int iii=0; iii < ncols*ncols; ++iii){
    if(sigma(iii) != sigmat(iii)){
      testsym = 1;
      break;
    }
  }
  
  if(testsym == 1){
    return arma::repmat(mu, 1, n) + (EigenvectorsSym(sigma) * arma::diagmat(sqrt(max(Eigenvalues(sigma), arma::zeros(ncols)))) * Y.t());
  }else{
    return arma::repmat(mu, 1, n) + (Eigenvectors(sigma) * arma::diagmat(sqrt(max(Eigenvalues(sigma), arma::zeros(ncols)))) * Y.t());
  }
}
//[[Rcpp::export]]
//Slice sampler for the latent lengths R in a projected normal regression model
arma::mat RSlice(arma::mat X1, arma::mat X2, arma::vec theta, arma::mat b1, arma::mat b2, int n, arma::mat r){
  
  arma::mat mub1 = b1*X1.t();
  arma::mat mub2 = b2*X2.t();
  arma::mat Dbd = cos(theta)%mub1.t() + sin(theta)%mub2.t();
  for (int jjj=0; jjj < n; ++jjj){
  arma::mat y = as<arma::vec>(Rcpp::runif(1,0,1)) % exp(-.5*pow((r.row(jjj)-Dbd.row(jjj)),2));
  arma::mat u = as<arma::vec>(Rcpp::runif(1,0,1));
  arma::mat r1 = Dbd.row(jjj) + max(-Dbd.row(jjj), -sqrt(-2*log(y)));
  arma::mat r2 = Dbd.row(jjj) + sqrt(-2*log(y));
  r.row(jjj)  = sqrt(((pow(r2,2)-pow(r1,2)) % u) + pow(r1,2));
}
  return r;
}

//[[Rcpp::export]]
//Sampling projected normal regression data with p predictors. Predictors are 
//centered after sampling outcome Y/Theta.
List RData(int N, int p, arma::vec Beta1, arma::vec Beta2, arma::vec mean, arma::vec sd){
  
  //Initialize matrices and vectors
  //Centered predictor matrices for both bivariate components.
  arma::mat X1(N, p);
  arma::mat X2(N, p);
  //Uncentered predictor matrices for both bivariate components.
  arma::mat X1old(N, p);
  arma::mat X2old(N, p);
  
  arma::mat Y(N,2);
  arma::vec Theta(N);
  arma::vec xc(N);
  arma::vec x(N);
  
  //Intercept
  arma::vec i = arma::ones<arma::vec>(N);
  
  X1old.col(0) = i;
  X2old.col(0) = i;
  X1.col(0) = i;
  X2.col(0) = i;
  
  //Sample predictors
  for(int kkk = 1; kkk < p; ++kkk){
    x = Rcpp::rnorm(N, mean(kkk-1), sd(kkk-1));
    X1old.col(kkk) = x;
    X2old.col(kkk) = x;
    X1.col(kkk) = x;
    X2.col(kkk) = x;
  }
  
  for(int iii = 0; iii < N; ++iii){
    
    xc(iii) = x(iii)-arma::mean(x);
    arma::vec mu(2);
    arma::mat x1 = X1.row(iii)*Beta1;
    arma::mat x2 = X2.row(iii)*Beta2;
    mu(0) = x1(0,0);
    mu(1) = x2(0,0);
    Y.row(iii) = mvrnormArmaEigen(1, mu, arma::eye(2,2)).t();
    double y1 = Y(iii,0);
    double y2 = Y(iii,1);
    Theta.row(iii) = atan2(y2,y1);
  }
  
  X1.col(1) = xc;
  X2.col(1) = xc;
  
  return  Rcpp::List::create(Rcpp::Named("Theta") = Theta,
                             Rcpp::Named("X1") = X1,
                             Rcpp::Named("X2") = X2,
			                       Rcpp::Named("X1old") = X1old,
                             Rcpp::Named("X2old") = X2old,
                             Rcpp::Named("Y") = Y);
}

// [[Rcpp::export]]
//Gibbs sampler for a projected normal regression model
List Regression(arma::vec theta, arma::mat X1, arma::mat X2, int tm, int tlag, int burn) {
  
  double n = theta.n_elem;
  double p1 = X1.n_cols;
  double p2 = X2.n_cols;
  
  arma::mat datose(n, 2);
  datose.col(0) = cos(theta);
  datose.col(1) = sin(theta);
  
  //Prior specification regression parameters
  arma::vec mu1 = arma::ones<arma::vec>(p1)*0;
  arma::vec mu2 = arma::ones<arma::vec>(p2)*0;
  arma::mat v1 = arma::eye<arma::mat>(p1,p1)*0.0001;
  arma::mat v2 = arma::eye<arma::mat>(p2,p2)*0.0001;
  
  //Posterior specification regression parameters
  arma::mat XtX1 = X1.t()*X1;
  arma::mat XtX2 = X2.t()*X2;
  arma::mat vstar1 = v1 + XtX1;
  arma::mat vstar2 = v2 + XtX2;
  arma::mat sigma1 = arma::inv(vstar1);
  arma::mat sigma2 = arma::inv(vstar2);
  arma::mat v1mu1 = v1*mu1;
  arma::mat v2mu2 = v2*mu2;
  
  //Initialize matrices for results
  int kk = tm*tlag;
  
  arma::mat r = arma::ones<arma::mat>(n,1);
  arma::mat B1(kk, p1);
  arma::mat B2(kk, p2);
  arma::mat Bc(kk, p1-1);
  arma::mat signeddistance(kk, p1-1);
  arma::mat ac(kk, p1-1);
  arma::mat ax(kk, p1-1);
  arma::mat mBc(kk, p1-1);
  arma::mat BcmX(kk, p1-1);
  arma::mat Y = r%datose.each_col();
  
  //burn-in
  for(int iii=0; iii < burn; ++iii){
 
  arma::mat XtY1 = X1.t()*Y.col(0); 
  arma::mat XtY2 = X2.t()*Y.col(1);
  arma::mat mstar1 = sigma1*(v1mu1 + XtY1);
  arma::mat mstar2 = sigma2*(v2mu2 + XtY2);
  
  //Sample coefficients
  arma::mat b1 = mvrnormArmaEigen(1, mstar1.col(0), sigma1).t();
  arma::mat b2 = mvrnormArmaEigen(1, mstar2.col(0), sigma2).t();
  
  //Sample R
  arma::mat mub1 = b1*X1.t();
  arma::mat mub2 = b2*X2.t();
  Rcpp::NumericMatrix Dbd = wrap(cos(theta)%mub1.t() + sin(theta)%mub2.t());
  r  = Dbd + (  pnorm(Dbd)/( dnorm(Dbd)+Dbd*pnorm(Dbd) )   );
  
  Y = r%datose.each_col();  
  
  }
  
  //Gibbs iterations
  for(int iii=0; iii < kk; ++iii){
    
  arma::mat XtY1 = X1.t()*Y.col(0);
  arma::mat XtY2 = X2.t()*Y.col(1);
  arma::mat mstar1 = sigma1*(v1mu1 + XtY1);
  arma::mat mstar2 = sigma2*(v2mu2 + XtY2);
  
  //Sample coefficients
  arma::mat b1 = mvrnormArmaEigen(1, mstar1.col(0), sigma1).t();
  arma::mat b2 = mvrnormArmaEigen(1, mstar2.col(0), sigma2).t();
  
  //Sample R
  r = RSlice(X1, X2, theta, b1, b2, n, r);
  
  Y = r%datose.each_col();
  
  //fill posterior coefficient matrices 
  B1.row(iii) = b1;
  B2.row(iii) = b2;
  for(int jjj=0; jjj < p1-1; ++jjj){
    ax.row(iii).col(jjj) = -((b1(0)*b1(jjj+1)+b2(0)*b2(jjj+1))/(pow(b1(jjj+1),2)+pow(b2(jjj+1),2)));
    double ac1 = b1(0)+b1(jjj+1)*ax.row(iii)(jjj);
    double ac2 = b2(0)+b2(jjj+1)*ax.row(iii)(jjj);
    ac.row(iii).col(jjj) =  atan2(ac2, ac1);
    signeddistance.row(iii).col(jjj) = sqrt(pow(b1(0)+b1(1)*ax.row(iii).col(jjj),2)+pow(b2(0)+b2(1)*ax.row(iii).col(jjj),2));
    double cutoff = atan2(b2(jjj+1), b1(jjj+1));
    double pi = arma::datum::pi;
    
    if(cutoff < 0){
      if((ac.row(iii)(jjj) < cutoff) | (ac.row(iii)(jjj) > (cutoff + pi))){
        signeddistance.row(iii).col(jjj) = -signeddistance.row(iii)(jjj);
      }else if((ac.row(iii)(jjj) > cutoff) & (ac.row(iii)(jjj) < (cutoff + pi))){
        signeddistance.row(iii).col(jjj) = signeddistance.row(iii)(jjj);
    }
    }
    
    if(cutoff > 0){
      if((ac.row(iii)(jjj) > cutoff) | (ac.row(iii)(jjj) < (cutoff - pi))){
        signeddistance.row(iii).col(jjj) = signeddistance.row(iii)(jjj);
      }else if((ac.row(iii)(jjj) < cutoff) & (ac.row(iii)(jjj) > (cutoff - pi))){
        signeddistance.row(iii).col(jjj) = -signeddistance.row(iii)(jjj);
    }
    }
    

    
    Bc.row(iii).col(jjj) =  tan(atan2(b2(0),b1(0)) - ac.row(iii)(jjj)) / -ax.row(iii)(jjj);
    BcmX.row(iii).col(jjj) =  Bc.row(iii)(jjj) / (1 + pow(Bc.row(iii)(jjj)*(mean(X1.col(jjj+1)) - ax.row(iii)(jjj)),2));
    arma::vec mBcP(n);
    for(int kkk=0; kkk < n; ++kkk){
      mBcP(kkk) = Bc.row(iii)(jjj) / (1 + pow(Bc.row(iii)(jjj)*(X1.row(kkk)(jjj+1) - ax.row(iii)(jjj)),2));
    }
    mBc.row(iii).col(jjj) =  mean(mBcP);
  }
 }
  
 arma::mat B = arma::join_rows(B1,B2);
   
   return  Rcpp::List::create(Rcpp::Named("B") = B,
                              Rcpp::Named("ax") = ax,
                              Rcpp::Named("ac") = ac,
                              Rcpp::Named("Bc") = Bc,
                              Rcpp::Named("BcmX") = BcmX,
                              Rcpp::Named("mBc") = mBc,
                              Rcpp::Named("SD") = signeddistance);
}

// [[Rcpp::export]]
//Projected normal regression simulation for 1 predictor
List SimRegRCPP1P(int N, int p , arma::vec Beta1, arma::vec Beta2, arma::vec mean, arma::vec sd, int tm, int tlag, int burn, int nsim, int EXBurn){
  
  //Initialize necessary matrixes
  arma::cube X1(N,p,nsim);
  arma::cube X2(N,p,nsim);
  arma::cube X1old(N,p,nsim);
  arma::cube X2old(N,p,nsim);
  arma::mat Theta(N,nsim);
  arma::cube Y(N,2,nsim);
  arma::cube B(tm*tlag, p*2, nsim);
  arma::cube Bc(tm*tlag, p-1, nsim);
  arma::cube ac(tm*tlag, p-1, nsim);
  arma::cube ax(tm*tlag, p-1, nsim);
  arma::cube signeddistance(tm*tlag, p-1, nsim);
  arma::cube mBc(tm*tlag, p-1, nsim);
  arma::cube BcmX(tm*tlag, p-1, nsim);
  
  arma::mat modeB(nsim,p*2);
  arma::vec modeax(nsim);
  arma::vec modeac(nsim);
  arma::vec modeBc(nsim);
  arma::vec modemBc(nsim);
  arma::vec modeBcmX(nsim);
  arma::vec modeSD(nsim);
  
  arma::vec InIntax(nsim);
  arma::vec InIntac(nsim);
  arma::vec InIntBc(nsim);
  arma::vec InIntmBc(nsim);
  arma::vec InIntBcmX(nsim);
  arma::vec ZInIntBc(nsim);
  arma::vec ZInIntmBc(nsim);
  arma::vec ZInIntBcmX(nsim);
  arma::mat InIntB(nsim, p*2);
  arma::vec InIntSD(nsim);
  arma::vec Accuracy(nsim);
  arma::vec NoEff(nsim);
  
  arma::vec biasax(nsim);
  arma::vec biasac(nsim);
  arma::vec biasBc(nsim);
  arma::vec biasmBc(nsim);
  arma::vec biasBcmX(nsim);
  arma::vec biasSD(nsim);

  arma::mat HPDax(nsim,2);
  arma::mat HPDac(nsim,2);
  arma::mat HPDBc(nsim,2);
  arma::mat HPDmBc(nsim,2);
  arma::mat HPDBcmX(nsim,2);
  arma::cube HPDB(nsim,2,p*2);
  arma::mat HPDSD(nsim,2);
  
  arma::vec realax(nsim);
  arma::vec realac(nsim);
  arma::vec realBc(nsim);
  arma::vec realmBc(nsim);
  arma::vec realBcmX(nsim);
  arma::vec realSD(nsim);
  
  //Initialize summary of results matrix
  Rcpp::NumericMatrix SumRes(p*2 + (p-1)*5 + 1, 6);
  Rcpp::List dimnms = 
  Rcpp::List::create(Rcpp::CharacterVector::create("B0I", "B1I", "B0II", "B1II", "ax", "ac", "Bc", "mBc", "BcmX", "SD"),
                       Rcpp::CharacterVector::create("RealValue", "EstMode", "Bias", "CILB", "CIHB", "Coverage"));
  SumRes.attr("dimnames") = dimnms;
  
  //Compute ax for population
  double axreal = -((Beta1(0)*Beta1(1)+Beta2(0)*Beta2(1))/(pow(Beta1(1),2)+pow(Beta2(1),2)));
  
  for(int iii=0; iii < nsim; ++iii){
    
    //Sample data
    Rcpp::Rcout << "Simulated Dataset: " << iii + 1 << std::endl;
    List Dataset = RData(N, p, Beta1, Beta2, mean, sd);
    X1.slice(iii) = Rcpp::as<arma::mat>(Dataset["X1"]);
    X2.slice(iii) = Rcpp::as<arma::mat>(Dataset["X2"]);
    X1old.slice(iii) = Rcpp::as<arma::mat>(Dataset["X1old"]);
    X2old.slice(iii) = Rcpp::as<arma::mat>(Dataset["X2old"]);
    Y.slice(iii) = Rcpp::as<arma::mat>(Dataset["Y"]);
    Theta.col(iii) = Rcpp::as<arma::vec>(Dataset["Theta"]);
    
    //Run the Gibbs sampler and save results
    List Results = Regression(Theta.col(iii), X1.slice(iii), X2.slice(iii), tm, tlag, burn);
  
    B.slice(iii) = Rcpp::as<arma::mat>(Results["B"]);
    Bc.slice(iii) = Rcpp::as<arma::mat>(Results["Bc"]);
    mBc.slice(iii) = Rcpp::as<arma::mat>(Results["mBc"]);
    BcmX.slice(iii) = Rcpp::as<arma::mat>(Results["BcmX"]);
    ac.slice(iii) = Rcpp::as<arma::mat>(Results["ac"]);
    ax.slice(iii) = Rcpp::as<arma::mat>(Results["ax"]);
    signeddistance.slice(iii) = Rcpp::as<arma::mat>(Results["SD"]);
    
    //Compute uncentered predictor mean
    double Xmean = arma::mean(X1old.slice(iii).col(1));
  
    //Compute real values
    //ax for individual datasets
    realax(iii) = -((Beta1(0)*Beta1(1)+Beta2(0)*Beta2(1))/(pow(Beta1(1),2)+pow(Beta2(1),2))) - Xmean;
    double realac1 = Beta1(0)+Beta1(1)*axreal;
    double realac2 = Beta2(0)+Beta2(1)*axreal;
    realac(iii) =  atan2(realac2, realac1);

    //Initialize Location/Accuracy indicator
    int Location = 1;
    
    //Compute more real values
    if(axreal == 0){
      realBc(iii) =  tan(atan2(Beta2(0)+Beta2(1),Beta1(0)+Beta1(1)) - realac(iii)) / (1-axreal);
    }
  
    if((Beta1(1)/Beta2(1)) - (Beta1(0)/Beta2(0)) == 0){
      realBc(iii) = 0;
      realBcmX(iii) = 0;
      realmBc(iii) = 0;
      Location = 0;
    }else if((Beta1(0) == Beta1(1)) & (Beta1(0) == 0)){
      realBc(iii) = 0;
      realBcmX(iii) = 0;
      realmBc(iii) = 0;
      Location = 0;
    }else if((Beta2(0) == Beta2(1)) & (Beta2(0) == 0)){
      realBc(iii) = 0;
      realBcmX(iii) = 0;
      realmBc(iii) = 0;
      Location = 0;
    }else if((Beta1(0)==Beta2(0)) & (Beta1(0) == 0)){
      realBc(iii) = 0;
      realBcmX(iii) = 0;
      realmBc(iii) = 0;
      Location = 0;
    }else if((Beta1(1)==Beta2(1)) & (Beta1(1) == 0)){
      realBc(iii) = 0;
      realBcmX(iii) = 0;
      realmBc(iii) = 0;
      Location = 0;
    }else if((Location == 1) & (axreal == 0)){
      realBc(iii) =  tan(atan2(Beta2(0)+Beta2(1),Beta1(0)+Beta1(1)) - realac(iii)) / (1-axreal);
      realBcmX(iii) =  realBc(iii) / (1 + pow(realBc(iii)*(arma::mean(X1.slice(iii).col(1)) - realax(iii)),2));
    
      arma::vec mBcP(N);
    
      for(int kkk=0; kkk < N; ++kkk){
        mBcP(kkk) = realBc(iii) / (1 + pow(realBc(iii)*(X1.slice(iii).row(kkk)(1) - realax(iii)),2));
      }
    
      realmBc(iii) =  arma::mean(mBcP);  
    }else{
      realBc(iii) =  tan(atan2(Beta2(0),Beta1(0)) - realac(iii)) / -axreal;
      realBcmX(iii) =  realBc(iii) / (1 + pow(realBc(iii)*(arma::mean(X1.slice(iii).col(1)) - realax(iii)),2));
    
      arma::vec mBcP(N);
    
      for(int kkk=0; kkk < N; ++kkk){
        mBcP(kkk) = realBc(iii) / (1 + pow(realBc(iii)*(X1.slice(iii).row(kkk)(1) - realax(iii)),2));
      }
    
      realmBc(iii) =  arma::mean(mBcP);
    }
  
  realSD(iii) = sqrt(pow(Beta1(0)+Beta1(1)*axreal,2)+pow(Beta2(0)+Beta2(1)*axreal,2));
    
    double cutoff = atan2(Beta2(1), Beta1(1));
    double pi = arma::datum::pi;
    
    if(cutoff < 0){
      if((realac(iii) < cutoff) | (realac(iii) > (cutoff + pi))){
        realSD(iii) = -realSD(iii);
      }else if((realac(iii) > cutoff) & (realac(iii) < (cutoff + pi))){
        realSD(iii) = realSD(iii);
      }
    }
    
    if(cutoff > 0){
      if((realac(iii) > cutoff) | (realac(iii) < (cutoff - pi))){
        realSD(iii) = realSD(iii);
      }else if((realac(iii) < cutoff) & (realac(iii) > (cutoff - pi))){
        realSD(iii) = -realSD(iii);
      }
    }
    
    //Compute summary statistics per dataset
    arma::vec MBc = Bc.slice(iii).col(0);
    arma::vec MmBc = mBc.slice(iii).col(0);
    arma::vec MBcmX = BcmX.slice(iii).col(0);
    arma::vec Max = ax.slice(iii).col(0);
    arma::vec Mac = ac.slice(iii).col(0);
    arma::vec Macmeas = signeddistance.slice(iii).col(0);
  
    modeBc(iii) = hmode(MBc.subvec(EXBurn-1, tm-1),0.1);
    modemBc(iii) = hmode(MmBc.subvec(EXBurn-1, tm-1),0.1);
    modeBcmX(iii) = hmode(MBcmX.subvec(EXBurn-1, tm-1),0.1);
    modeax(iii) = hmode(Max.subvec(EXBurn-1, tm-1), 0.1);
    modeac(iii)= hmodeC(Mac.subvec(EXBurn-1, tm-1), 0.1);
    modeSD(iii) = hmode(Macmeas.subvec(EXBurn-1, tm-1), 0.1);
  
    biasax(iii) = realax(iii)-modeax(iii);
    biasac(iii) = realac(iii)-modeac(iii);
    biasBc(iii) = realBc(iii)-modeBc(iii);
    biasmBc(iii) = realmBc(iii)-modemBc(iii);
    biasBcmX(iii) = realBcmX(iii)-modeBcmX(iii);
    biasSD(iii) = realSD(iii)-modeSD(iii);

    HPDBc.row(iii) = hmodeci(MBc.subvec(EXBurn-1, tm-1),0.95).t();
    HPDmBc.row(iii) = hmodeci(MmBc.subvec(EXBurn-1, tm-1),0.95).t();
    HPDBcmX.row(iii) = hmodeci(MBcmX.subvec(EXBurn-1, tm-1),0.95).t();
    HPDax.row(iii) = hmodeci(Max.subvec(EXBurn-1, tm-1), 0.95).t();
    HPDac.row(iii) = hmodeciC(Mac.subvec(EXBurn-1, tm-1), 0.95).t();
    HPDSD.row(iii) = hmodeci(Macmeas.subvec(EXBurn-1, tm-1), 0.95).t();
    
    //Check whether real values and 0 are in HPD interval
    if((HPDBc.row(iii)(0) < realBc(iii)) & (HPDBc.row(iii)(1) > realBc(iii))){
      InIntBc(iii) = 1;
    }else{
      InIntBc(iii) = 0;
    }
  
    if((HPDmBc.row(iii)(0) < realmBc(iii)) & (HPDmBc.row(iii)(1) > realmBc(iii))){
      InIntmBc(iii) = 1;
    }else{
      InIntmBc(iii) =0;
    }
  
    if((HPDBcmX.row(iii)(0) < realBcmX(iii)) & (HPDBcmX.row(iii)(1) > realBcmX(iii))){
      InIntBcmX(iii) = 1;
    }else{
      InIntBcmX(iii) = 0;
    }
  
    if((HPDax.row(iii)(0) < realax(iii)) & (HPDax.row(iii)(1) > realax(iii))){
      InIntax(iii) = 1;
    }else{
      InIntax(iii) = 0;
    }
    
    if((HPDSD.row(iii)(0) < realSD(iii)) & (HPDSD.row(iii)(1) > realSD(iii))){
      InIntSD(iii) = 1;
    }else{
      InIntSD(iii) =0;
    }
    
    if((HPDSD.row(iii)(0) < 0) & (HPDSD.row(iii)(1) > 0)){
      Accuracy(iii) = 1;
    }else{
      Accuracy(iii) = 0;
    }
  
    if((HPDBc.row(iii)(0) < 0) & (HPDBc.row(iii)(1) > 0)){
      ZInIntBc(iii) = 1;
    }else{
      ZInIntBc(iii) = 0;
    }
  
    if((HPDmBc.row(iii)(0) < 0) & (HPDmBc.row(iii)(1) > 0)){
      ZInIntmBc(iii) = 1;
    }else{
      ZInIntmBc(iii) = 0;
    }
  
    if((HPDBcmX.row(iii)(0) < 0) & (HPDBcmX.row(iii)(1) > 0)){
      ZInIntBcmX(iii) = 1;
    }else{
      ZInIntBcmX(iii) =0;
    }
  
    if(HPDac.row(iii)(0) > HPDac.row(iii)(1)){
        if((realac(iii) > HPDac.row(iii)(0)) | ((0 <= realac(iii)) & (realac(iii) < HPDac.row(iii)(1)))){
          InIntac(iii) = 1;
        }else{
          InIntac(iii) = 0;
        }
      
      }else{
        
          if((realac(iii) > HPDac.row(iii)(0)) & (realac(iii) < HPDac.row(iii)(1))){
            InIntac(iii) = 1;
          }else{
            InIntac(iii) = 0;
          }
      }
    
    for(int jjj=0; jjj < p*2; ++jjj){
    
      //Real values
      SumRes(0,0) = Beta1(0);
      SumRes(1,0) = Beta1(1);
      SumRes(2,0) = Beta2(0);
      SumRes(3,0) = Beta2(1);
    
      arma::vec MB = B.slice(iii).col(jjj);
    
      modeB(iii,jjj) = hmode(MB.subvec(EXBurn-1, tm-1), 0.1);
    
      HPDB.slice(jjj).row(iii) = hmodeci(MB.subvec(EXBurn-1, tm-1), 0.95).t();
      
      
      if((HPDB.slice(jjj).row(iii)(0) < SumRes(jjj,0)) & (HPDB.slice(jjj).row(iii)(1) > SumRes(jjj,0))){
        InIntB(iii,jjj) = 1;
        }else{
        InIntB(iii,jjj) = 0;
        }
      
    }
    
    if((HPDB.slice(1).row(iii)(0) < 0) & (HPDB.slice(1).row(iii)(1) > 0) & (HPDB.slice(3).row(iii)(0) < 0) & (HPDB.slice(3).row(iii)(1) > 0)){
      NoEff(iii) = 1;
    }else{
      NoEff(iii) = 0;
    }
  
  }
  
  //Compute summary statistics for simulation design per parameter
  for(int jjj=0; jjj < p*2; ++jjj){
    SumRes(jjj,1) = arma::mean(modeB.col(jjj));
    SumRes(jjj,2) = SumRes(jjj,0)-SumRes(jjj,1);
    SumRes(jjj,3) = arma::mean(HPDB.slice(jjj).col(0));
    SumRes(jjj,4) = arma::mean(HPDB.slice(jjj).col(1));
    SumRes(jjj,5) = sum(InIntB.col(jjj))/nsim;
  }
  
  SumRes(p*2, 1) = arma::mean(modeax);
  SumRes(p*2, 0) = arma::mean(realax);   
  SumRes(p*2, 2) = arma::mean(biasax);
  SumRes(p*2, 3) = arma::mean(HPDax.col(0));
  SumRes(p*2, 4) = arma::mean(HPDax.col(1));
  SumRes(p*2, 5) = sum(InIntax)/nsim;
  
  SumRes(p*2 + 1, 1) = hmodeC(modeac, 0.1);
  SumRes(p*2 + 1, 0) = arma::mean(realac);
  SumRes(p*2 + 1, 2) = arma::mean(biasac);
  SumRes(p*2 + 1, 3) = hmodeC(HPDac.col(0), 0.1);
  SumRes(p*2 + 1, 4) = hmodeC(HPDac.col(1), 0.1);
  SumRes(p*2 + 1, 5) = sum(InIntac)/nsim;
  
  SumRes(p*2 + 2, 1) = arma::mean(modeBc);
  SumRes(p*2 + 2, 0) = arma::mean(realBc);
  SumRes(p*2 + 2, 2) = arma::mean(biasBc);
  SumRes(p*2 + 2, 3) = arma::mean(HPDBc.col(0));
  SumRes(p*2 + 2, 4) = arma::mean(HPDBc.col(1));
  SumRes(p*2 + 2, 5) = sum(InIntBc)/nsim;
  
  SumRes(p*2 + 3, 1) = arma::mean(modemBc);
  SumRes(p*2 + 3, 0) = arma::mean(realmBc);
  SumRes(p*2 + 3, 2) = arma::mean(biasmBc);
  SumRes(p*2 + 3, 3) = arma::mean(HPDmBc.col(0));
  SumRes(p*2 + 3, 4) = arma::mean(HPDmBc.col(1));
  SumRes(p*2 + 3, 5) = sum(InIntmBc)/nsim;
  
  SumRes(p*2 + 4, 1) = arma::mean(modeBcmX);
  SumRes(p*2 + 4, 0) = arma::mean(realBcmX);
  SumRes(p*2 + 4, 2) = arma::mean(biasBcmX);
  SumRes(p*2 + 4, 3) = arma::mean(HPDBcmX.col(0));
  SumRes(p*2 + 4, 4) = arma::mean(HPDBcmX.col(1));
  SumRes(p*2 + 4, 5) = sum(InIntBcmX)/nsim;
  
  SumRes(p*2 + 5, 1) = arma::mean(modeSD);
  SumRes(p*2 + 5, 0) = arma::mean(realSD);
  SumRes(p*2 + 5, 2) = arma::mean(biasSD);
  SumRes(p*2 + 5, 3) = hmode(HPDSD.col(0), 0.1);
  SumRes(p*2 + 5, 4) = hmode(HPDSD.col(1), 0.1);
  SumRes(p*2 + 5, 5) = sum(InIntSD)/nsim;
  
  return  Rcpp::List::create(Rcpp::Named("SumRes") = SumRes,
                             Rcpp::Named("Theta") = Theta,
                             Rcpp::Named("X1") = X1,
                             Rcpp::Named("X2") = X2,
                             Rcpp::Named("X1old") = X1old,
                             Rcpp::Named("X2old") = X2old,
                             Rcpp::Named("Y") = Y,
                             Rcpp::Named("B") = B,
                             Rcpp::Named("ax") = ax,
                             Rcpp::Named("ac") = ac,
                             Rcpp::Named("SD") = signeddistance,
                             Rcpp::Named("Bc") = Bc,
                             Rcpp::Named("BcmX") = BcmX,
                             Rcpp::Named("mBc") = mBc,
                             Rcpp::Named("NoEffindicator") = NoEff,
                             Rcpp::Named("Accuracyindicator") = Accuracy,
                             Rcpp::Named("ZeroEffBc") = ZInIntBc,
                             Rcpp::Named("ZeroEffmBc") = ZInIntmBc,
                             Rcpp::Named("ZeroEffBcmX") = ZInIntBcmX);
}