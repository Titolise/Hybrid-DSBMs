// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// Likelihood
// [[Rcpp::export]]

List complk_dyn_cpp(arma::cube Y, arma::vec piv, arma::mat Pi, arma::mat Psi, arma::umat Ub, int k, int n, int TT) {
  arma::cube Phi(n, k, TT);
  for(int t=0; t<TT; t++) {
    arma::uvec Ub_t = Ub.col(t);
    arma::mat Y_t = Y.slice(t);
    for(int i = 0; i < n; i++) {
      for(int u = 0; u < k; u++) {
        double log_phi = 0.0;
        for(int j = 0; j < n; j++) {
          double y_ij = Y_t(i, j);
          if(i == j || std::isnan(y_ij)) continue;
          int v = Ub_t(j) - 1;
          double psi_uv = std::max(std::min(Psi(u, v), 1.0 - 1e-16), 1e-16);
          log_phi += y_ij * std::log(psi_uv) + (1.0 - y_ij) * std::log(1.0 - psi_uv);
        }
        Phi(i, u, t) = std::max(std::exp(log_phi), 1e-300);
      }
    }
  }

  arma::cube L(n, k, TT, arma::fill::zeros);
  for(int i=0; i<n; i++) {
    for(int u=0; u<k; u++) L(i, u, 0) = Phi(i, u, 0) * piv(u);
  }
  for(int t=1; t<TT; t++) {
    L.slice(t) = Phi.slice(t) % (L.slice(t-1) * Pi);
  }

  arma::vec pv = sum(L.slice(TT-1), 1);
  double lk = 0;
  for(int i=0; i<n; i++) lk += std::log(std::max(pv(i), 1e-300));

  return List::create(Named("lk") = lk, Named("Phi") = Phi, Named("L") = L, Named("pv") = pv);
}

// Greedy Search Originale in C++
// [[Rcpp::export]]

List greedy_search_cpp_original(arma::cube Y, arma::umat Ub, arma::mat Psi, arma::mat Pi, arma::vec piv) {
  int n = Y.n_rows;
  int TT = Y.n_slices;
  int k = Psi.n_rows;

  arma::cube Phi(n, k, TT);
  for(int t=0; t<TT; t++) {
    arma::uvec Ub_t = Ub.col(t);
    arma::mat Y_t = Y.slice(t);
    for(int i = 0; i < n; i++) {
      for(int u = 0; u < k; u++) {
        double log_phi = 0.0;
        for(int j = 0; j < n; j++) {
          double y_ij = Y_t(i, j);
          if(i == j || std::isnan(y_ij)) continue;
          int v = Ub_t(j) - 1;
          double psi_uv = std::max(std::min(Psi(u, v), 1.0 - 1e-16), 1e-16);
          log_phi += y_ij * std::log(psi_uv) + (1.0 - y_ij) * std::log(1.0 - psi_uv);
        }
        Phi(i, u, t) = std::max(std::exp(log_phi), 1e-300);
      }
    }
  }

  arma::cube L(n, k, TT);
  for(int i=0; i<n; i++) {
    for(int u=0; u<k; u++) L(i, u, 0) = Phi(i, u, 0) * piv(u);
  }
  for(int t=1; t<TT; t++) L.slice(t) = Phi.slice(t) % (L.slice(t-1) * Pi);

  arma::vec pv = sum(L.slice(TT-1), 1);
  double best_global_lk = 0;
  for(int i=0; i<n; i++) best_global_lk += std::log(std::max(pv(i), 1e-300));

  for(int t=0; t<TT; t++) {
    for(int h=0; h<n; h++) {
      int current_class = Ub(h, t) - 1;
      int best_class = current_class;
      double best_lk = best_global_lk;
      arma::mat best_Phi_t = Phi.slice(t);

      for(int v=0; v<k; v++) {
        if(v == current_class) continue;
        Ub(h, t) = v + 1;

        arma::mat Phi_temp_t = Phi.slice(t);
        arma::uvec Ub_t = Ub.col(t);
        arma::mat Y_t = Y.slice(t);

        for(int i = 0; i < n; i++) {
          for(int u = 0; u < k; u++) {
            double log_phi = 0.0;
            for(int j = 0; j < n; j++) {
              double y_ij = Y_t(i, j);
              if(i == j || std::isnan(y_ij)) continue;
              int c_v = Ub_t(j) - 1;
              double psi_uv = std::max(std::min(Psi(u, c_v), 1.0 - 1e-16), 1e-16);
              log_phi += y_ij * std::log(psi_uv) + (1.0 - y_ij) * std::log(1.0 - psi_uv);
            }
            Phi_temp_t(i, u) = std::max(std::exp(log_phi), 1e-300);
          }
        }

        arma::cube L_temp = L;
        if(t == 0) {
          for(int i=0; i<n; i++) {
            for(int u=0; u<k; u++) L_temp(i, u, 0) = Phi_temp_t(i, u) * piv(u);
          }
          for(int tt=1; tt<TT; tt++) L_temp.slice(tt) = Phi.slice(tt) % (L_temp.slice(tt-1) * Pi);
        } else {
          L_temp.slice(t) = Phi_temp_t % (L_temp.slice(t-1) * Pi);
          for(int tt=t+1; tt<TT; tt++) L_temp.slice(tt) = Phi.slice(tt) % (L_temp.slice(tt-1) * Pi);
        }

        arma::vec pv_temp = sum(L_temp.slice(TT-1), 1);
        double lk_temp = 0;
        for(int i=0; i<n; i++) lk_temp += std::log(std::max(pv_temp(i), 1e-300));

        if(std::isfinite(lk_temp) && lk_temp > best_lk) {
          best_lk = lk_temp;
          best_class = v;
          best_Phi_t = Phi_temp_t;
        }
        Ub(h, t) = current_class + 1;
      }

      if(best_class != current_class) {
        Ub(h, t) = best_class + 1;
        best_global_lk = best_lk;
        Phi.slice(t) = best_Phi_t;

        if(t == 0) {
          for(int i=0; i<n; i++) {
            for(int u=0; u<k; u++) L(i, u, 0) = Phi(i, u, 0) * piv(u);
          }
          for(int tt=1; tt<TT; tt++) L.slice(tt) = Phi.slice(tt) % (L.slice(tt-1) * Pi);
        } else {
          L.slice(t) = Phi.slice(t) % (L.slice(t-1) * Pi);
          for(int tt=t+1; tt<TT; tt++) L.slice(tt) = Phi.slice(tt) % (L.slice(tt-1) * Pi);
        }
      }
    }
  }
  pv = sum(L.slice(TT-1), 1);
  return List::create(Named("lk") = best_global_lk, Named("Phi") = Phi, Named("L") = L, Named("pv") = pv, Named("Ub") = Ub);
}
