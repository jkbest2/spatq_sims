#include <TMB.hpp>

//' Log-normal density
//'
//' Calculate the density of an observation that is log-normally distributed
//' with a given mean and standard deviation of the distribution on the log
//' scale. Vectorized to accept vectors<Type> observations, means, and/or
//' standard deviations.
//' @param x [Type] observation
//' @param meanlog [Type] mean on the log scale
//' @param sdlog [Type] standard deviation on the log scale
//' @param give_log [int] return density (0, default) or log density (1)
//' @return Type density or log density
template <class Type>
Type dlnorm(Type x, Type meanlog, Type sdlog, int give_log = 0) {
  Type logres = dnorm(log(x), meanlog, sdlog, true) - log(x);
  if (give_log)
    return logres;
  else
    return exp(logres);
}
VECTORIZE4_ttti(dlnorm);

template <class Type> Type objective_function<Type>::operator()() {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // DATA section
  // Vector of observed catches; zero or positive
  DATA_VECTOR(obs);        // 2,000
  // Abundance projection matrices
  DATA_SPARSE_MATRIX(A);   // 2,000 × 404
  DATA_SPARSE_MATRIX(IA);  // 10,000 × 404
  // FEM matrices for SPDE spatial and spatiotemporal effects. Sharing the mesh
  // between effects means that this only needs to be passed once.
  DATA_STRUCT(spde, spde_t);
  DATA_INTEGER(link_fn);

  // PARAMETER section
  // Abundance fixed effects
  PARAMETER_VECTOR(beta1);        // n_year
  PARAMETER_VECTOR(beta2);        // n_year
  // Abundance spatial effects
  PARAMETER_VECTOR(omega1);       // 404
  PARAMETER_VECTOR(omega2);       // 404
  // Spatial field parameters
  PARAMETER_VECTOR(log_kappa);    // 2
  PARAMETER_VECTOR(log_tau);      // 2
  // Log catch variation parameter
  PARAMETER(log_sigma);

  //  Negative log-likelihood
  vector<Type> nll(3);
  nll.setZero();
  // Number of observations
  int n_obs = obs.size();
  // Number of index stations
  int n_ind = IA.rows();
  // Spatial field parameters
  vector<Type> kappa = exp(log_kappa);
  vector<Type> itau = exp(-log_tau);
  Type sigma = exp(log_sigma);

  // Spatial field contribution to the likelihood
  SparseMatrix<Type> Q1 = Q_spde(spde, kappa(0));
  nll(0) += GMRF(Q1)(omega1);
  SparseMatrix<Type> Q2 = Q_spde(spde, kappa(1));
  nll(1) += GMRF(Q2)(omega2);

  //  Spatial field
  vector<Type> spat1(n_obs);
  spat1 = A * omega1;
  vector<Type> spat2(n_obs);
  spat2 = A * omega2;

  // Expected values
  vector<Type> p1(n_obs);
  p1 = beta1(0) + itau(0) * spat1;
  vector<Type> p2(n_obs);
  p2 = beta2(0) + itau(1) * spat2;

  vector<Type> p_enc(n_obs);
  vector<Type> pos_r(n_obs);
  for (int i = 0; i < n_obs; i++) {
    if (link_fn == 0) {
      p_enc(i) = invlogit(p1(i));
      pos_r(i) = exp(p2(i));
    } else if (link_fn == 1) {
      p_enc(i) = Type(1.0) - exp(-exp(p1(i)));
      pos_r(i) = exp(p1(i) + p2(i)) / p_enc(i);
    }

    if (obs(i) == 0) {
      nll(2) -= log(1 - p_enc(i));
    } else {
      nll(2) -= log(p_enc(i));
      nll(2) -= dlnorm(obs(i), log(pos_r(i)) - sigma * sigma / Type(2.0), sigma, true);
    }
  }

  REPORT(p_enc);
  REPORT(pos_r);
  REPORT(nll);

  return(nll.sum());
}
