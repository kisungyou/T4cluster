#include "RcppArmadillo.h"
#include "utilities.h"
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Repair 1: stable sphere primitives; no float conversion or small-step cutoff.
namespace {

const double sphere_antipodal_tol = 1e-12;

double sphere_norm(const arma::rowvec& x) {
  if (!x.is_finite()) {
    Rcpp::stop("* sphere : a vector contains nonfinite values.");
  }
  const double scale = arma::abs(x).max();
  if (scale == 0.0) return 0.0;
  return scale * arma::norm(x / scale, 2);
}

arma::rowvec sphere_unit(const arma::rowvec& x) {
  if (x.n_elem < 2 || !x.is_finite()) {
    Rcpp::stop("* sphere : a direction must have at least two finite coordinates.");
  }
  const double scale = arma::abs(x).max();
  if (scale == 0.0) {
    Rcpp::stop("* sphere : a zero vector has no direction.");
  }
  arma::rowvec scaled = x / scale;
  return scaled / arma::norm(scaled, 2);
}

struct SpherePair {
  double angle;
  double axial;
  double tangent_norm;
  double radius;
  arma::rowvec tangent;
};

SpherePair sphere_pair(const arma::rowvec& x, const arma::rowvec& y) {
  if (x.n_elem != y.n_elem || x.n_elem < 2 ||
      !x.is_finite() || !y.is_finite()) {
    Rcpp::stop("* sphere : incompatible or nonfinite directions.");
  }
  const double xx = arma::dot(x, x);
  const double xy = arma::dot(x, y);
  if (!(xx > 0.0) || !std::isfinite(xx) || !std::isfinite(xy)) {
    Rcpp::stop("* sphere : invalid base direction.");
  }
  SpherePair pair;
  pair.axial = xy / std::sqrt(xx);
  pair.tangent = y - (xy / xx) * x;
  pair.tangent -= (arma::dot(x, pair.tangent) / xx) * x;
  pair.tangent_norm = sphere_norm(pair.tangent);
  pair.radius = std::hypot(pair.axial, pair.tangent_norm);
  if (!(pair.radius > 0.0) || !std::isfinite(pair.radius)) {
    Rcpp::stop("* sphere : invalid target direction.");
  }
  pair.angle = std::atan2(pair.tangent_norm, pair.axial);
  return pair;
}

bool sphere_antipodal(const SpherePair& pair) {
  return pair.axial < 0.0 &&
    pair.tangent_norm <= sphere_antipodal_tol * pair.radius;
}

arma::rowvec sphere_pair_log(const SpherePair& pair) {
  if (sphere_antipodal(pair)) {
    Rcpp::stop("* gskmeans : a required logarithm is at a numerical antipode; choose another initialization or partition.");
  }
  if (pair.tangent_norm == 0.0) {
    return arma::rowvec(pair.tangent.n_elem, arma::fill::zeros);
  }
  // Divide the direction before multiplying by the angle to avoid overflow.
  return pair.angle * (pair.tangent / pair.tangent_norm);
}

} // anonymous namespace

arma::rowvec geometry_proj(arma::rowvec x, arma::rowvec z) {
  const double xx = arma::dot(x, x);
  if (!(xx > 0.0) || !std::isfinite(xx)) {
    Rcpp::stop("* sphere : invalid projection base.");
  }
  return z - (arma::dot(x, z) / xx) * x;
}

double geometry_dist(arma::rowvec x, arma::rowvec y) {
  return sphere_pair(x, y).angle;
}

arma::rowvec geometry_exp(arma::rowvec x, arma::rowvec u) {
  x = sphere_unit(x);
  if (u.n_elem != x.n_elem || !u.is_finite()) {
    Rcpp::stop("* sphere : invalid tangent vector.");
  }
  u = geometry_proj(x, u);
  const double radius = sphere_norm(u);
  if (radius == 0.0) return x;
  if (!std::isfinite(radius)) {
    Rcpp::stop("* sphere : nonfinite tangent norm.");
  }
  double sinc;
  if (radius < 1e-4) {
    const double r2 = radius * radius;
    sinc = 1.0 - r2 / 6.0 + (r2 * r2) / 120.0;
  } else {
    sinc = std::sin(radius) / radius;
  }
  return sphere_unit(std::cos(radius) * x + sinc * u);
}

arma::rowvec geometry_log(arma::rowvec x, arma::rowvec y) {
  return sphere_pair_log(sphere_pair(x, y));
}

// Repair 2: shared input, label, and fixed-K empty-cluster checks.
namespace {

arma::mat sphere_data(const arma::mat& X, int K, int maxiter,
                      double abstol, const std::string& initializer) {
  if (X.n_rows == 0 || X.n_cols < 2 || !X.is_finite()) {
    Rcpp::stop("* spherical k-means : data must be a nonempty finite matrix with at least two columns.");
  }
  if (X.n_rows > static_cast<arma::uword>(std::numeric_limits<int>::max())) {
    Rcpp::stop("* spherical k-means : too many rows for the current initialization helpers.");
  }
  if (K < 1 || static_cast<arma::uword>(K) > X.n_rows) {
    Rcpp::stop("* spherical k-means : k must lie between one and the number of observations.");
  }
  if (maxiter < 1 || !std::isfinite(abstol) || !(abstol > 0.0)) {
    Rcpp::stop("* spherical k-means : maxiter must be positive and abstol must be finite and positive.");
  }
  if (initializer != "kmeans" && initializer != "gmm") {
    Rcpp::stop("* spherical k-means : initializer must be 'kmeans' or 'gmm'.");
  }
  // Normalize a local copy, including when a native caller bypasses the R wrapper.
  arma::mat data = X;
  for (arma::uword i = 0; i < data.n_rows; ++i) {
    data.row(i) = sphere_unit(data.row(i));
  }
  return data;
}

arma::field<arma::uvec> sphere_indices(const arma::uvec& labels,
                                    arma::uword N, int K,
                                    const std::string& stage) {
  if (labels.n_elem != N || labels.n_elem == 0 ||
      labels.max() >= static_cast<arma::uword>(K)) {
    Rcpp::stop("* spherical k-means : invalid labels during %s.", stage.c_str());
  }
  arma::field<arma::uvec> indices(K);
  for (int k = 0; k < K; ++k) {
    indices(k) = arma::find(labels == static_cast<arma::uword>(k));
    if (indices(k).n_elem == 0) {
      Rcpp::stop("* spherical k-means : empty cluster during %s; try another initialization or a smaller k.", stage.c_str());
    }
  }
  return indices;
}

arma::uvec sphere_initialize(const arma::mat& X, int K,
                            const std::string& initializer, int maxiter) {
  if (K == 1) return arma::uvec(X.n_rows, arma::fill::zeros);
  if (initializer == "kmeans") {
    return arma::trans(label_kmeans(X, K, maxiter));
  }
  // The requested initializer is called exactly once.
  return arma::trans(label_gmm(X, K, maxiter));
}

bool sphere_same_labels(const arma::uvec& a, const arma::uvec& b) {
  return a.n_elem == b.n_elem && arma::accu(a != b) == 0;
}

double sphere_objective_slack(arma::uword N, double a, double b) {
  return 128.0 * std::numeric_limits<double>::epsilon() *
    (static_cast<double>(N) + std::abs(a) + std::abs(b) + 1.0);
}

} // anonymous namespace

// Repair 3: maximize dot products for both spherical objectives.
arma::uvec sp_spkmeans_label(arma::mat X, arma::mat centers) {
  if (centers.n_rows == 0 || centers.n_cols != X.n_cols ||
      !centers.is_finite()) {
    Rcpp::stop("* spherical k-means : invalid centers in assignment.");
  }
  arma::uvec labels(X.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < X.n_rows; ++i) {
    double best = -std::numeric_limits<double>::infinity();
    arma::uword winner = 0;
    for (arma::uword k = 0; k < centers.n_rows; ++k) {
      const double score = arma::dot(X.row(i), centers.row(k));
      if (score > best) { // Exact ties retain the first center.
        best = score;
        winner = k;
      }
    }
    labels(i) = winner;
  }
  return labels;
}

arma::uvec sp_gskmeans_label(arma::mat X, arma::mat centers) {
  return sp_spkmeans_label(X, centers);
}

// Repair 4: exact zero-resultant convention; positive resultants stay exact.
arma::mat sp_spkmeans_centers(arma::mat X, arma::field<arma::uvec> index) {
  arma::mat centers(index.n_elem, X.n_cols, arma::fill::zeros);
  for (arma::uword k = 0; k < index.n_elem; ++k) {
    const arma::uvec& rows = index(k);
    if (rows.n_elem == 0) {
      Rcpp::stop("* spkmeans : empty cluster in center update.");
    }
    arma::rowvec resultant = arma::sum(X.rows(rows), 0);
    if (!resultant.is_finite()) {
      Rcpp::stop("* spkmeans : nonfinite cluster resultant.");
    }
    if (arma::abs(resultant).max() == 0.0) {
      // Any unit center minimizes the fixed-partition loss in this case.
      centers.row(k) = sphere_unit(X.row(rows(0)));
    } else {
      centers.row(k) = sphere_unit(resultant);
    }
  }
  return centers;
}

// Both costs evaluate the supplied partition AND supplied centers. This matters
// after final reassignment if the iteration budget ended before a fixed point.
double sp_spkmeans_cost(arma::mat X, arma::mat centers,
                       arma::field<arma::uvec> index) {
  double score = 0.0;
  for (arma::uword k = 0; k < index.n_elem; ++k) {
    const arma::uvec& rows = index(k);
    for (arma::uword i = 0; i < rows.n_elem; ++i) {
      score += arma::dot(X.row(rows(i)), centers.row(k));
    }
  }
  return score; // Higher is better; equals summed resultant lengths at centroids.
}

double sp_gskmeans_cost(arma::mat X, arma::mat centers,
                       arma::field<arma::uvec> index) {
  double cost = 0.0;
  for (arma::uword k = 0; k < index.n_elem; ++k) {
    const arma::uvec& rows = index(k);
    for (arma::uword i = 0; i < rows.n_elem; ++i) {
      const double angle = geometry_dist(X.row(rows(i)), centers.row(k));
      cost += angle * angle;
    }
  }
  return cost; // Lower is better; unhalved sum of squared geodesic distances.
}

// Repair 5: warm-started intrinsic means with residual-based stopping.
// A small average-log residual certifies stationarity, not a global minimum.
namespace {

struct SphereMeanResult {
  arma::rowvec center;
  bool converged;
  bool stalled;
  int iterations;
  double residual;
};

bool sphere_mean_evaluate(const arma::mat& X, const arma::rowvec& center,
                          double& objective, arma::rowvec& direction) {
  objective = 0.0;
  direction.zeros(X.n_cols);
  for (arma::uword i = 0; i < X.n_rows; ++i) {
    const SpherePair pair = sphere_pair(center, X.row(i));
    // Trials at numerical antipodes can be rejected without throwing.
    if (sphere_antipodal(pair)) return false;
    objective += 0.5 * pair.angle * pair.angle;
    direction += sphere_pair_log(pair);
  }
  objective /= static_cast<double>(X.n_rows);
  direction /= static_cast<double>(X.n_rows);
  direction = geometry_proj(center, direction);
  return std::isfinite(objective) && direction.is_finite();
}

SphereMeanResult sphere_frechet(const arma::mat& X,
                               const arma::rowvec& initial,
                               int mean_maxiter, double mean_tol) {
  SphereMeanResult result;
  result.center = sphere_unit(initial);
  result.converged = false;
  result.stalled = false;
  result.iterations = 0;
  double value;
  arma::rowvec direction(X.n_cols, arma::fill::zeros);
  if (!sphere_mean_evaluate(X, result.center, value, direction)) {
    Rcpp::stop("* gskmeans : a required logarithm is at a numerical antipode; choose another initialization or partition.");
  }
  result.residual = sphere_norm(direction);
  result.converged = result.residual <= mean_tol;
  for (int iteration = 0; iteration < mean_maxiter && !result.converged; ++iteration) {
    const double squared_residual = result.residual * result.residual;
    const double slack = 64.0 * std::numeric_limits<double>::epsilon() *
      (1.0 + std::abs(value));
    bool accepted = false;
    double beta = 0.5;
    for (int trial = 0; trial < 25; ++trial, beta *= 0.5) {
      arma::rowvec candidate = geometry_exp(result.center, beta * direction);
      double candidate_value;
      arma::rowvec candidate_direction(X.n_cols, arma::fill::zeros);
      if (!sphere_mean_evaluate(X, candidate, candidate_value, candidate_direction)) {
        continue;
      }
      // Slack softens Armijo only; the evaluated objective may NOT increase.
      if (candidate_value <= value &&
          candidate_value <= value - 1e-4 * beta * squared_residual + slack) {
        const double movement = sphere_norm(candidate - result.center);
        result.center = candidate;
        value = candidate_value;
        direction = candidate_direction;
        result.residual = sphere_norm(direction);
        result.iterations += 1;
        result.converged = result.residual <= mean_tol;
        result.stalled = movement == 0.0 && !result.converged;
        accepted = true;
        break;
      }
    }
    if (!accepted) result.stalled = true;
    if (result.stalled) break;
  }
  return result;
}

arma::mat sphere_intrinsic_centers(const arma::mat& X,
                                   const arma::field<arma::uvec>& index,
                                   const arma::mat& initial,
                                   int mean_maxiter, double mean_tol,
                                   Rcpp::IntegerVector& iterations,
                                   Rcpp::LogicalVector& converged,
                                   Rcpp::LogicalVector& stalled,
                                   Rcpp::NumericVector& residual) {
  arma::mat centers = initial;
  for (arma::uword k = 0; k < index.n_elem; ++k) {
    Rcpp::checkUserInterrupt();
    if (index(k).n_elem == 0) {
      Rcpp::stop("* gskmeans : empty cluster in center update.");
    }
    const SphereMeanResult result = sphere_frechet(
      X.rows(index(k)), initial.row(k), mean_maxiter, mean_tol);
    centers.row(k) = result.center;
    iterations[k] = result.iterations;
    converged[k] = result.converged;
    stalled[k] = result.stalled;
    residual[k] = result.residual;
  }
  return centers;
}

bool sphere_all_true(const Rcpp::LogicalVector& x) {
  for (R_xlen_t i = 0; i < x.size(); ++i) if (x[i] != TRUE) return false;
  return true;
}

bool sphere_any_true(const Rcpp::LogicalVector& x) {
  for (R_xlen_t i = 0; i < x.size(); ++i) if (x[i] == TRUE) return true;
  return false;
}

// Repair 6: one shared outer loop and consistent final labels/cost/diagnostics.
Rcpp::List sphere_fit(const arma::mat& X, int K, const std::string& initializer,
                      int maxiter, double abstol, bool printer,
                      bool intrinsic, int mean_maxiter, double mean_tol) {
  arma::mat data = sphere_data(X, K, maxiter, abstol, initializer);
  if (intrinsic && (mean_maxiter < 1 || !std::isfinite(mean_tol) || !(mean_tol > 0.0))) {
    Rcpp::stop("* gskmeans : mean_maxiter must be positive and mean_tol must be finite and positive.");
  }
  const std::string name = intrinsic ? "gskmeans" : "spkmeans";
  arma::uvec labels = sphere_initialize(data, K, initializer, maxiter);
  arma::field<arma::uvec> index = sphere_indices(labels, data.n_rows, K, "initialization");
  arma::mat centers = sp_spkmeans_centers(data, index);
  Rcpp::IntegerVector mean_iterations(K, 0);
  Rcpp::LogicalVector mean_converged(K, true);
  Rcpp::LogicalVector mean_stalled(K, false);
  Rcpp::NumericVector mean_residual(K, 0.0);
  if (intrinsic) {
    centers = sphere_intrinsic_centers(data, index, centers, mean_maxiter,
      mean_tol, mean_iterations, mean_converged, mean_stalled, mean_residual);
  }
  double cost = intrinsic ? sp_gskmeans_cost(data, centers, index) :
    sp_spkmeans_cost(data, centers, index);
  if (!std::isfinite(cost)) Rcpp::stop("* spherical k-means : nonfinite initial objective.");
  std::vector<double> history;
  history.push_back(cost);
  bool converged = false;
  bool stopped_for_mean = false;
  int iterations = 0;

  for (int iteration = 0; iteration < maxiter; ++iteration) {
    Rcpp::checkUserInterrupt();
    arma::uvec candidate_labels = sp_spkmeans_label(data, centers);
    arma::field<arma::uvec> candidate_index = sphere_indices(
      candidate_labels, data.n_rows, K, "reassignment");
    const bool same_partition = sphere_same_labels(candidate_labels, labels);
    arma::mat candidate_centers;
    if (intrinsic) {
      candidate_centers = sphere_intrinsic_centers(data, candidate_index, centers,
        mean_maxiter, mean_tol, mean_iterations, mean_converged,
        mean_stalled, mean_residual);
    } else {
      candidate_centers = sp_spkmeans_centers(data, candidate_index);
    }
    const double candidate_cost = intrinsic ?
      sp_gskmeans_cost(data, candidate_centers, candidate_index) :
      sp_spkmeans_cost(data, candidate_centers, candidate_index);
    if (!std::isfinite(candidate_cost) || !candidate_centers.is_finite()) {
      Rcpp::stop("* spherical k-means : nonfinite center or objective during iteration.");
    }
    const double slack = sphere_objective_slack(data.n_rows, cost, candidate_cost);
    const double improvement = intrinsic ? cost - candidate_cost : candidate_cost - cost;
    if (improvement < -slack) {
      Rcpp::stop("* spherical k-means : an update worsened the objective beyond numerical tolerance.");
    }
    const double change = std::abs(candidate_cost - cost);
    labels = candidate_labels;
    index = candidate_index;
    centers = candidate_centers;
    cost = candidate_cost;
    iterations = iteration + 1;
    history.push_back(cost);
    if (printer) {
      Rcpp::Rcout << "* " << name << " : iteration " << iterations
                  << " : objective " << cost << std::endl;
    }

    if (same_partition && change <= abstol &&
        (!intrinsic || sphere_all_true(mean_converged))) {
      // A small intrinsic refinement can move a Voronoi boundary. Verify the
      // candidate centers themselves before reporting a fixed partition.
      const arma::uvec check_labels = sp_spkmeans_label(data, centers);
      sphere_indices(check_labels, data.n_rows, K, "convergence check");
      if (sphere_same_labels(check_labels, labels)) {
        converged = true;
        break;
      }
    }
    if (intrinsic && same_partition && sphere_any_true(mean_stalled)) {
      stopped_for_mean = true;
      break;
    }
  }

  // Final reassignment does not perform another center update. If labels
  // change, means need not be centroids of the returned partition.
  const arma::uvec final_labels = sp_spkmeans_label(data, centers);
  const bool labels_stable = sphere_same_labels(final_labels, labels);
  index = sphere_indices(final_labels, data.n_rows, K, "final reassignment");
  const double final_cost = intrinsic ? sp_gskmeans_cost(data, centers, index) :
    sp_spkmeans_cost(data, centers, index);
  if (!std::isfinite(final_cost)) Rcpp::stop("* spherical k-means : nonfinite final objective.");
  const double final_improvement = intrinsic ? cost - final_cost : final_cost - cost;
  if (final_improvement < -sphere_objective_slack(data.n_rows, cost, final_cost)) {
    Rcpp::stop("* spherical k-means : final reassignment worsened the objective beyond numerical tolerance.");
  }
  cost = final_cost;
  history.push_back(cost); // Always append, even if unchanged: length=iterations+2.

  if (intrinsic) {
    // Residuals and convergence flags describe the RETURNED partition.
    // mean_iterations still records the last centroid solve before this refresh.
    for (int k = 0; k < K; ++k) {
      double value;
      arma::rowvec direction(data.n_cols, arma::fill::zeros);
      if (!sphere_mean_evaluate(data.rows(index(k)), centers.row(k), value, direction)) {
        Rcpp::stop("* gskmeans : a logarithm required for final residuals is at a numerical antipode.");
      }
      mean_residual[k] = sphere_norm(direction);
      mean_converged[k] = mean_residual[k] <= mean_tol;
    }
  }
  converged = converged && labels_stable &&
    (!intrinsic || sphere_all_true(mean_converged));
  std::string status = "converged";
  if (!converged) {
    status = stopped_for_mean ? "mean_stalled" : "maxiter";
    Rcpp::warning("* %s : convergence conditions were not met (status=%s). Returned labels and cost were recomputed at the returned means.",
                  name.c_str(), status.c_str());
  }
  if (printer && converged) {
    Rcpp::Rcout << "* " << name << " : convergence achieved." << std::endl;
  }

  Rcpp::List output = Rcpp::List::create(
    Rcpp::Named("means") = centers,
    Rcpp::Named("cluster") = final_labels,
    Rcpp::Named("cost") = cost,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("status") = status,
    Rcpp::Named("labels_stable") = labels_stable,
    Rcpp::Named("objective_history") = Rcpp::wrap(history));
  if (intrinsic) {
    output["mean_converged"] = mean_converged;
    output["mean_iterations"] = mean_iterations;
    output["mean_residual"] = mean_residual;
  }
  return output;
}

} // anonymous namespace

// Repair 7: preserve the existing six-argument native entry points.
// [[Rcpp::export]]
Rcpp::List sp_spkmeans(arma::mat& X, int K, std::string initializer,
                      int maxiter, double abstol, bool printer) {
  return sphere_fit(X, K, initializer, maxiter, abstol, printer, false, 100, 1e-8);
}

Rcpp::List sp_gskmeans_impl(arma::mat& X, int K, std::string initializer,
                          int maxiter, double abstol, bool printer,
                          int mean_maxiter, double mean_tol) {
  return sphere_fit(X, K, initializer, maxiter, abstol, printer, true,
                    mean_maxiter, mean_tol);
}

// [[Rcpp::export]]
Rcpp::List sp_gskmeans(arma::mat& X, int K, std::string initializer,
                      int maxiter, double abstol, bool printer) {
  return sp_gskmeans_impl(X, K, initializer, maxiter, abstol, printer, 100, 1e-8);
}

// [[Rcpp::export]]
Rcpp::List sp_gskmeans_control(arma::mat& X, int K, std::string initializer,
                              int maxiter, double abstol, bool printer,
                              int mean_maxiter = 100, double mean_tol = 1e-8) {
  return sp_gskmeans_impl(X, K, initializer, maxiter, abstol, printer,
                         mean_maxiter, mean_tol);
}
