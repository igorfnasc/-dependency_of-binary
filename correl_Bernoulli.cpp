#include <cmath>
#include <Rmath.h>
#include "RcppArmadillo.h"
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

using namespace arma;
using namespace Rcpp;

#include <RcppArmadilloExtensions/sample.h>

typedef std::vector<double> stdvec;


// [[Rcpp::export]]
// arma::mat
arma::mat func_bernoulli_correl_fam(arma::mat matrixprobs,
                                              int ks,
                                              double rho1,
                                              double rho2,
                                              int kk){
  
  
  // List ret;
  int tamj = matrixprobs.n_rows;
  arma::mat vectabuasimulada(kk,tamj);

  vectabuasimulada.fill(0);
  
  for(int j = 0 ; j < (kk);j++){
    for(int i = 0 ; i < (ks);i++){
      int    vivo1 = 1 ;
      int    vivo2 = 1 ;
      int    vivo3 = 1 ;
      
      int contador = 0;
      int alea1=0;
      int alea2=0;
      int alea3=0;
      
      double p1 = 0;
      double p2 = 0;
      double p3 = 0;
      
      
      
      while((vivo1+vivo2+vivo3)>0){
        
        double p1t = matrixprobs(contador,0);
        double p2t = matrixprobs(contador,1);      
        double p3t = matrixprobs(contador,2);      
        
        if(vivo1==0){
          p1t = 1.0;
        }
        
        if(vivo2==0){
          p2t = 1.0;
        }
        
        if(vivo3==0){
          p3t = 1.0;
        }
        
        
          alea1 = R::rbinom( 1,p1t);
          
          if(alea1==1){
            vivo1=0;
          }
          
          // # conjuge
          double p2tt = p2t + rho1*((alea1)*(1-p2t) - (1-alea1)*(p2t));
          
          
          if(vivo2==0){
            p2tt = 1.0;
          }
          
          
          alea2 = R::rbinom( 1,p2tt);
          // alea2 = 1;
          
          if(alea2==1){
            vivo2=0;
          }
          
          
        // # filho
        double p3tt = p3t + rho2*( (alea1)*(1-p3t) - (1-alea1)*(p3t)  );
          
          if(vivo3==0){
            p3tt = 1.0;
          }
          
          
          alea3 = R::rbinom( 1,p3tt);
          // alea2 = 1;
          
          if(alea3==1){
            vivo3=0;
          }
        

        
        vectabuasimulada(j,contador)= vectabuasimulada(j,contador)+  max(NumericVector::create(vivo1,vivo2,vivo3));
        contador = contador + 1;
        
      }
    }
    
  }
  
  return(vectabuasimulada);
}



// [[Rcpp::export]]
// arma::mat
arma::mat func_bernoulli_correl_coup(arma::mat matrixprobs,
                                              int ks,
                                              double rho1,
                                              int kk){
  
  
  // List ret;
  int tamj = matrixprobs.n_rows;
  arma::mat vectabuasimulada(kk,tamj);
  arma::mat ps(2,1);
  arma::mat vs(2,1);
  
  vectabuasimulada.fill(0);
  
  for(int j = 0 ; j < (kk);j++){
    for(int i = 0 ; i < (ks);i++){
      int    vivo1 = 1 ;
      int    vivo2 = 1 ;
      
      int   contador = 0;
      int alea1=0;
      int alea2=0;
      
      double p1 = 0;
      double p2 = 0;
      
      
      
      while(vivo1==1|vivo2==1){
        
        double p1t = matrixprobs(contador,0);
        double p2t = matrixprobs(contador,1);      
        
        if(vivo1==0){
          p1t = 1.0;
        }
        
        if(vivo2==0){
          p2t = 1.0;
        }
        
        
        alea1 = R::rbinom( 1,p1t);
        
        if(alea1==1){
          vivo1=0;
        }
        
        // # mulher
        double p2tt = p2t + rho1*( (alea1)*(1-p2t) - (1-alea1)*(p2t)  );
        
        
        
        alea2 = R::rbinom( 1,p2tt);
        // alea2 = 1;
        
        if(alea2==1){
          vivo2=0;
        }
        
        
        
        
        vectabuasimulada(j,contador)= vectabuasimulada(j,contador)+ max(NumericVector::create(vivo1,vivo2));
        // vectabuasimulada(i,contador)=  alea1;
        // vectabuasimulada(0,contador)= alea2;
        contador = contador + 1;
        
      }
    }
    
  }
  
  return(vectabuasimulada);
}



