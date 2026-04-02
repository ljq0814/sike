#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <cmath>
#define _USE_MATH_DEFINES
#include <algorithm>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace std;
using namespace Rcpp;
#define INFINITYNUMBER 1E7

// [[Rcpp::export]]
arma::vec spsign(arma::vec x){
    unsigned p = x.n_elem;
    double nx = arma::norm(x,2);
    if (nx < 1e-3){
        arma::vec res(p,arma::fill::zeros);
        return(res);
    }
    arma::vec res = x/nx;
    return(res);
}

// [[Rcpp::export]]
arma::mat spcov(arma::mat X){
    unsigned n = X.n_rows;unsigned p = X.n_cols;
    arma::mat res(p,p,arma::fill::zeros);
    for(unsigned i = 0; i < n; ++i){
        arma::rowvec tmp = X.row(i);
        double nx = arma::norm(tmp,2);
        if (nx < 1e-3){
            continue;
        }
        res = res + tmp.t()*(tmp)/pow(nx,2)/n;
    }
    return(res);
}

// [[Rcpp::export]]
arma::mat kendalltau(arma::mat X){
    unsigned n = X.n_rows;unsigned p = X.n_cols;
    arma::mat res(p,p,arma::fill::zeros);
    for(unsigned i = 0; i < n; ++i){
        for (unsigned j = i+1; j < n; ++j){
            arma::rowvec tmp = X.row(i) - X.row(j);
            double nx = arma::norm(tmp,2);
            if (nx < 1e-3){
                continue;
            }
            res = res + tmp.t()*(tmp)/pow(nx,2)/n/(n-1)*2;
        }
    }
    return(res);
}

arma::vec updateloc(arma::mat resi, arma::mat sqB){
    unsigned n = resi.n_rows;unsigned p = resi.n_cols;
    arma::vec res(p,p,arma::fill::zeros);
    arma::vec nume(p, arma::fill::zeros);
    double deno = 0.0;
    for (unsigned i = 0; i < n; ++i){
        nume = nume + spsign(arma::conv_to<arma::vec>::from(resi.row(i)));
        deno = deno + 1/arma::norm(resi.row(i),"2");
    }
    return(sqB*nume/deno);
}

arma::vec updateshape(arma::mat resi, arma::mat sqB){
    unsigned n = resi.n_rows;unsigned p = resi.n_cols;
    arma::vec res(p,p,arma::fill::zeros);
    for (unsigned i = 0; i < n; ++i){
        arma::vec ue = spsign(arma::conv_to<arma::vec>::from(resi.row(i)));
        res = res + ue*ue.t();
    }
    return(sqB*res*sqB/n*p);
}

// [[Rcpp::export]]
arma::mat hrstandardization(arma::mat X, double epsilon = 5e-3, unsigned maxiter = 1000){
    unsigned n = X.n_rows;unsigned p = X.n_cols;
    arma::mat initcov = kendalltau(X)/p;
    arma::vec initloc(p,arma::fill::zeros);
    double criteria = 1E7;
    unsigned i = 0;
    arma::mat resi(p,p,arma::fill::zeros);

    while (criteria > epsilon && i < maxiter){
        arma::vec eigcovval;
        arma::mat eigcovvec;
	    eig_sym(eigcovval, eigcovvec, initcov);
	    
        arma::mat sqB = eigcovvec*(arma::diagmat(arma::sqrt(arma::clamp(eigcovval,0.0,1E7))))*(eigcovvec.t());
        arma::mat invsqB = eigcovvec*(arma::diagmat(arma::pow(arma::clamp(eigcovval,0.0,1E7),-0.5)))*(eigcovvec.t());
        
        arma::vec allone(n, arma::fill::ones);
        resi = (X - allone*initloc.t())*invsqB;
        initloc = updateloc(resi,sqB);
        initcov = updateshape(resi,sqB);

    }
    return(initcov);
}



