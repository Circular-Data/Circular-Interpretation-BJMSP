// [[Rcpp::depends(RcppArmadillo)]]

#include <math.h>
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::vec Eigenvalues(arma::mat X) {
  
  arma::cx_vec eigval;
  arma::eig_gen(eigval, X);
  
  return arma::conv_to<arma::vec>::from(eigval);
}

// [[Rcpp::export]]
arma::mat Eigenvectors(arma::mat X) {
  
  arma::cx_vec eigval;
  arma::cx_mat eigvec;
  arma::eig_gen(eigval, eigvec, X);
  
  arma::mat eigvc = arma::conv_to<arma::mat>::from(eigvec);
  return  -1 * eigvc;
}

// [[Rcpp::export]]

arma::mat mvrnormArmaEigen(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y(n, ncols);
  for(int iii=0; iii < ncols; ++iii){
    Y.col(iii) = as<arma::vec>(rnorm(n));
  }
  return arma::repmat(mu, 1, n) +(Eigenvectors(sigma) * arma::diagmat(sqrt(max(Eigenvalues(sigma), arma::zeros(ncols)))) * Y.t());
}

//[[Rcpp::export]]

arma::mat RSlice(arma::mat X1, arma::mat X2, arma::vec theta, arma::mat b1, arma::mat b2, int n, arma::mat r){
  
  arma::mat mub1 = b1*X1.t();
  arma::mat mub2 = b2*X2.t();
  arma::mat Dbd = cos(theta)%mub1.t() + sin(theta)%mub2.t();
  for (int jjj=0; jjj < n; ++jjj){
    arma::mat y = as<arma::vec>(runif(1,0,1)) % exp(-.5*pow((r.row(jjj)-Dbd.row(jjj)),2));
    arma::mat u = as<arma::vec>(runif(1,0,1));
    arma::mat r1 = Dbd.row(jjj) + max(-Dbd.row(jjj), -sqrt(-2*log(y)));
    arma::mat r2 = Dbd.row(jjj) + sqrt(-2*log(y));
    r.row(jjj)  = sqrt(((pow(r2,2)-pow(r1,2)) % u) + pow(r1,2));
  }
  return r;
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