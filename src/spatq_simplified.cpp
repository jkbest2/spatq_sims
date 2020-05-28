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

//' log-Poisson link
//'
//' Calculate the log-probability of zero cathch, encounter log-probability and
//' log-positive catch rate using the Poisson link introduced in Thorson (2017).
//' @param log_n [Type] log group density
//' @param log_w [Type] log weight per group
//' @param a area sampled offset, defaults to 1.0
//' @return A vector of length 3 containing the log-probability of encounter,
//'   log-probability of zero catch, and log-positive catch rate (particularly
//'   useful when using a log-normal catch distribution).
template <class Type>
vector<Type> logpoislink(Type log_n, Type log_w, Type a = Type(1.0)) {
  // Initialize 3-vector that will contain probability of zero catch,
  // probability of encounter, and positive catch rate
  // TODO Can this be done in-place by passing in a view on a matrix or similar?
  vector<Type> log_ppr(3);

  // First calculate log-probability of zero catch, as it is used to calculated
  // the probability of zero catch.
  log_ppr(1) = -a * exp(log_n);
  // First argument is `log(1) == 0`
  log_ppr(0) = logspace_sub(Type(0.0), log_ppr(1));
  log_ppr(2) = log_n - log_ppr(0) + log_w;

  return log_ppr;
}

//' Random log-normal deviates
//'
//' Generate independent log-normal draws with mean and standard deviation on
//' the log scale. Length of output is the max of length of meanlog or sdlog.
//' @param meanlog mean on the log scale
//' @param meansd standard deviation on the log scale
//' @return log-normally distributed random deviate
// template<class Type>
// Type rlnorm(Type meanlog, Type sdlog) {
//   return exp(rnorm(meanlog, sdlog));
// }
// VECTORIZE2_tt(rlnorm);
// VECTORIZE2_n(rlnorm);
template <class Type> Type rlnorm(Type meanlog, Type sdlog) {
  return exp(Rf_rnorm(asDouble(meanlog), asDouble(sdlog)));
}
VECTORIZE2_tt(rlnorm);
VECTORIZE2_n(rlnorm);

template <class Type> Type objective_function<Type>::operator()() {
  using namespace R_inla;
  using namespace density;
  using namespace Eigen;

  // ===========================================================================
  // DATA section
  // ---------------------------------------------------------------------------
  // Vector of observed catches; zero or positive
  DATA_VECTOR(catch_obs);
  // Area swept offset for each observation
  DATA_VECTOR(area_swept);

  // Abundance fixed effects design matrix
  DATA_MATRIX(X_n);
  DATA_MATRIX(X_w);
  // Index fixed effects design matrix
  DATA_MATRIX(IX_n);
  DATA_MATRIX(IX_w);

  // Abundance projection matrices
  DATA_SPARSE_MATRIX(A_spat);   // Spatial
  DATA_SPARSE_MATRIX(A_sptemp); // Spatiotemporal
  // Index projection matrices
  DATA_SPARSE_MATRIX(IA_spat);   // Spatial
  DATA_SPARSE_MATRIX(IA_sptemp); // Spatiotemporal

  // Index integration weights
  DATA_VECTOR(Ih);

  // ---------------------------------------------------------------------------
  // Fixed effects design matrix
  DATA_MATRIX(R_n);
  DATA_MATRIX(R_w);

  // ---------------------------------------------------------------------------
  // FEM matrices for SPDE spatial and spatiotemporal effects. Sharing the mesh
  // between effects means that this only needs to be passed once.
  DATA_STRUCT(spde, spde_t);

  // ---------------------------------------------------------------------------
  // Vector indicating which of the 12 random processes should be included.
  // Currently takes a length-6 vector, and can only switch off pairs of numbers
  // density and weight-per-group processes. Indexing is currently:
  // 0, 1: omega_n, omega_w
  // 2, 3: epsilon1_n, epsilon1_w
  DATA_IVECTOR(proc_switch);

  // ---------------------------------------------------------------------------
  // Flags to control GMRF normalization and return early for normalization
  DATA_INTEGER(norm_flag);
  DATA_INTEGER(incl_data);

  // ===========================================================================
  // PARAMETER section
  // ---------------------------------------------------------------------------
  // Abundance fixed effects
  PARAMETER_VECTOR(beta_n); // number of fixed effects
  PARAMETER_VECTOR(beta_w); // number of fixed effects

  // Abundance spatial effects
  PARAMETER_VECTOR(omega_n); // N_vert
  PARAMETER_VECTOR(omega_w); // N_vert

  // Abundance spatiotemporal effects
  // PARAMETER_MATRIX(epsilon1_n); // N_vert × N_yrs - 1
  // PARAMETER_MATRIX(epsilon1_w); // N_vert × N_yrs - 1
  PARAMETER_MATRIX(epsilon_n); // N_vert × N_yrs
  PARAMETER_MATRIX(epsilon_w); // N_vert × N_yrs

  // --------------------------------------------------------------------------
  // Catchability fixed effects
  PARAMETER_VECTOR(lambda_n); // number of fixed effects
  PARAMETER_VECTOR(lambda_w); // number of fixed effects

  // Spatial and spatiotemporal field parameters
  PARAMETER_VECTOR(log_kappa); // 2
  PARAMETER_VECTOR(log_tau);   // 2

  // Log catch variation parameter
  PARAMETER(log_sigma); // 1

  // ===========================================================================
  // Derived values
  // ---------------------------------------------------------------------------
  // Get number of observations
  int N_obs = catch_obs.size();
  // Get number of years
  int N_yrs = epsilon_n.cols();
  // Get number of integration locations for each index year
  int N_I = Ih.size();
  // Convert norm_flag incl_data to boolean
  bool nrmlz = bool(norm_flag);
  bool incl_datalik = bool(incl_data);

  // Put unconstrained parameters on their natural (constrained) scales
  vector<Type> kappa = exp(log_kappa);
  vector<Type> itau = exp(-log_tau);
  Type sigma = exp(log_sigma);

  // ===========================================================================
  // Log-likelihood accumulator
  // ---------------------------------------------------------------------------
  // Initialize joint negative log-likelihood accumulator; accumulate spatial
  // and spatiotemporal effects in indices corresponding to indices of kappa and
  // tau.
  // 0,1: Abundance spatial
  // 2,3: Abundance spatiotemporal
  // 4:  Observation likelihood
  vector<Type> jnll(5);
  jnll.setZero();

  // ===========================================================================
  // Abundance fixed effects
  // ---------------------------------------------------------------------------
  vector<Type> fixef_n(N_obs);
  vector<Type> fixef_w(N_obs);
  fixef_n = X_n * beta_n;
  fixef_w = X_w * beta_w;

  // Index fixed effects
  vector<Type> Ifixef_n(N_I);
  vector<Type> Ifixef_w(N_I);
  Ifixef_n = IX_n * beta_n;
  Ifixef_w = IX_w * beta_w;

  // ===========================================================================
  // Abundance spatial effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  vector<Type> spat_n(N_obs);
  vector<Type> spat_w(N_obs);
  vector<Type> Ispat_n(N_I);
  vector<Type> Ispat_w(N_I);
  // Initialize them to zeros
  spat_n.setZero();
  spat_w.setZero();
  Ispat_n.setZero();
  Ispat_w.setZero();

  if (proc_switch(0)) {
    spat_n = A_spat * omega_n;
    SparseMatrix<Type> Q_n_om = Q_spde(spde, kappa(0));
    SCALE_t<GMRF_t<Type>> gmrf_n_om = SCALE(GMRF(Q_n_om, nrmlz), itau(0));
    jnll(0) += gmrf_n_om(omega_n);

    Ispat_n = IA_spat * omega_n;
  }
  if (proc_switch(1)) {
    spat_w = A_spat * omega_w;
    SparseMatrix<Type> Q_w_om = Q_spde(spde, kappa(1));
    SCALE_t<GMRF_t<Type>> gmrf_w_om = SCALE(GMRF(Q_w_om, nrmlz), itau(1));
    jnll(1) += gmrf_w_om(omega_w);

    Ispat_w = IA_spat * omega_w;
  }

  // ===========================================================================
  // Abundance spatiotemporal effects
  // ---------------------------------------------------------------------------
  // Project spatial effects from mesh nodes to observation locations
  // matrix<Type> epsilon_n(epsilon1_n.rows(), N_yrs);
  // matrix<Type> epsilon_w(epsilon1_w.rows(), N_yrs);
  vector<Type> sptemp_n(N_obs);
  vector<Type> sptemp_w(N_obs);
  vector<Type> Isptemp_n(N_I);
  vector<Type> Isptemp_w(N_I);
  // epsilon_n.setZero();
  // epsilon_w.setZero();
  sptemp_n.setZero();
  sptemp_w.setZero();
  Isptemp_n.setZero();
  Isptemp_w.setZero();

  if (proc_switch(2)) {
    // epsilon_n.leftCols(N_yrs - 1) = epsilon1_n;
    // epsilon_n.rightCols(1) = -epsilon1_n.rowwise().sum();

    sptemp_n = A_sptemp * epsilon_n.value();
    sptemp_n *= itau(2);

    SparseMatrix<Type> Q_n_ep = Q_spde(spde, kappa(2));
    GMRF_t<Type> gmrf_n_ep = GMRF(Q_n_ep, nrmlz);

    for (int yr = 0; yr < N_yrs; yr++) {
      jnll(2) += gmrf_n_ep(epsilon_n.col(yr));
    }

    REPORT(epsilon_n);

    Isptemp_n = IA_sptemp * epsilon_n.value();
    Isptemp_n *= itau(2);
  }
  if (proc_switch(3)) {
    // epsilon_w.leftCols(N_yrs - 1) = epsilon1_w;
    // epsilon_w.rightCols(1) = -epsilon1_w.rowwise().sum();

    sptemp_w = A_sptemp * epsilon_w.value();
    sptemp_w *= itau(3);

    SparseMatrix<Type> Q_w_ep = Q_spde(spde, kappa(3));
    GMRF_t<Type> gmrf_w_ep = GMRF(Q_w_ep, nrmlz);

    for (int yr = 0; yr < N_yrs; yr++) {
      jnll(3) += gmrf_w_ep(epsilon_w.col(yr));
    }

    REPORT(epsilon_w);

    Isptemp_w = IA_sptemp * epsilon_w.value();
    Isptemp_w *= itau(3);
  }

  // ===========================================================================
  // Return early if not including data likelihood
  // ---------------------------------------------------------------------------
  if (!incl_datalik)
    return (jnll.sum());

  // ===========================================================================
  // Catchability fixed effects
  // ---------------------------------------------------------------------------
  vector<Type> qfixef_n(N_obs);
  vector<Type> qfixef_w(N_obs);
  qfixef_n = R_n * lambda_n;
  qfixef_w = R_w * lambda_w;

  // ===========================================================================
  // Calculate linear predictor
  // ---------------------------------------------------------------------------
  // Get group density (n) and weight per group (w) for each observation
  vector<Type> log_n(N_obs);
  vector<Type> log_w(N_obs);
  log_n = fixef_n + spat_n + sptemp_n + qfixef_n;
  log_w = fixef_w + spat_w + sptemp_w + qfixef_w;

  // Index linear predictor
  vector<Type> Ilog_n(N_I);
  vector<Type> Ilog_w(N_I);
  Ilog_n = Ifixef_n + Ispat_n + Isptemp_n;
  Ilog_w = Ifixef_w + Ispat_w + Isptemp_w;

  // ===========================================================================
  // Apply link function
  // ---------------------------------------------------------------------------
  // Get group density (n) and weight per group (w) for each observation
  // Convert to log-probability encounter (p), log-probability of zero catch
  // (1mp) and positive catch rate (r)
  // TODO Use a N_obs×3 matrix here instead of three vectors. Maybe transpose
  // that so that elements to be accessed together are in a column?
  vector<Type> log_p_enc(N_obs);
  vector<Type> log_p_zero(N_obs);
  vector<Type> log_r(N_obs);
  // Temporary vector to use for return of logpoislink function. See todo above
  // for probably more elegant solution.
  vector<Type> log_ppr_tmp(3);
  for (int i = 0; i < N_obs; i++) {
    // TODO VECTORIZE2_tt poislink function? Not sure if this will work if it
    // returns a vector. Probably better to pass in views to p and r here.
    log_ppr_tmp = logpoislink(log_n(i), log_w(i), area_swept(i));
    log_p_enc(i) = log_ppr_tmp(0);
    log_p_zero(i) = log_ppr_tmp(1);
    log_r(i) = log_ppr_tmp(2);
  }

  // ===========================================================================
  // Observation likelihood
  // ---------------------------------------------------------------------------
  for (int i = 0; i < N_obs; i++) {
    if (catch_obs(i) == 0) {
      jnll(4) -= log_p_zero(i);
    } else {
      jnll(4) -= log_p_enc(i) + dlnorm(catch_obs(i),
                                       log_r(i) - sigma * sigma / Type(2.0),
                                       sigma, true);
    }
  }

  // ===========================================================================
  // Index of abundance
  // ---------------------------------------------------------------------------
  vector<Type> Index(N_yrs);
  Index.setZero();
  int i0 = 0;
  int i1 = 0;

  int N_I_per_yr = N_I / N_yrs;

  for (int yr = 0; yr < N_yrs; yr++) {
    i0 = yr * N_I_per_yr;
    i1 = (yr + 1) * N_I_per_yr;
    for (int i = i0; i < i1; i++) {
      // Should prevent some numerical issues; cribbed from VAST index
      // calculation.
      Index(yr) += Ih(i) * exp(Ilog_n(i)) * exp(Ilog_w(i));
    }
  }

  // ===========================================================================
  // Reports
  // ---------------------------------------------------------------------------
  vector<Type> rho_sp;
  vector<Type> sigma_sp;
  rho_sp = sqrt(8) / kappa;
  sigma_sp = itau / (kappa * 2 * sqrt(PI));

  REPORT(jnll);
  REPORT(Ilog_n);
  REPORT(Ilog_w);
  REPORT(rho_sp);
  REPORT(sigma_sp);

  ADREPORT(Index);

  return jnll.sum();
}
