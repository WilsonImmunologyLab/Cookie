#include <Rcpp.h>
using namespace Rcpp;


//' binary (XOR) distance coding for char and bool factors
//'
//' @param data data frame
//' @export
// [[Rcpp::export]]
NumericMatrix binaryCodingCpp(DataFrame data) {
  int n = data.nrows();
  int m = data.length();
  NumericMatrix res(n,n);
  double cutoff = 10;
  Rcout << "Start running... \n";
  for (int i = 0; i < n - 1; i++) {
    double progress = 100*i/n;
    if(progress > cutoff) {
      Rcout << cutoff << "% processed... \n";
      cutoff = cutoff + 10;
    }

    for (int j = i + 1; j < n; j++) {
      int vec_res = 0;
      for (int k = 0; k < m; k++) {
        CharacterVector working_vec = data[k];
        if(working_vec[i] != "NA" && working_vec[j] != "NA") {
          if(working_vec[i] != working_vec[j]) {
            vec_res++;
          }
        }
      }
      res(i,j) = vec_res;
      res(j,i) = vec_res;
    }
  }
  return res;
}

//' hamming distance coding for numerical factors
//'
//' @param data data frame
//' @export
// [[Rcpp::export]]
NumericMatrix hammingCodingCpp(DataFrame data) {
  int n = data.nrows();
  int m = data.length();
  NumericMatrix res(n,n);
  double cutoff = 10;
  Rcout << "Start running... \n";
  for (int i = 0; i < n - 1; i++) {
    double progress = 100*i/n;
    if(progress > cutoff) {
      Rcout << cutoff << "% processed... \n";
      cutoff = cutoff + 10;
    }

    for (int j = i + 1; j < n; j++) {
      int vec_res = 0;
      for (int k = 0; k < m; k++) {
        NumericVector working_vec = data[k];
        if(working_vec[i] != -1 && working_vec[j] != -1) {
          vec_res = vec_res + abs(working_vec[i] - working_vec[j]);
        }
      }
      res(i,j) = vec_res;
      res(j,i) = vec_res;
    }
  }
  return res;
}
