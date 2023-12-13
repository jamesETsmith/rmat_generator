#ifndef ____RMAT_UTILS_HPP
#define ____RMAT_UTILS_HPP

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>  // std::setprecision
#include <iostream>
#include <mutex>
#include <numeric>
#include <thread>
#include <vector>

void print_matrix(std::vector<double>& mat) {
  size_t n = static_cast<size_t>(std::lround(std::sqrt(mat.size())));
  assert(n * n == mat.size());
  n = std::min(n, 16UL);

  //   std::cout;
  for (size_t i = 0; i < n; i++) {
    for (size_t j = 0; j < n; j++) {
      std::cout << std::fixed << mat[i * n + j] << " ";
    }
    std::cout << std::endl;
  }
}

void kronecker_prod(std::vector<double>& c, std::vector<double>& a,
                    std::vector<double>& b) {
  // assert(a.size() == b.size());
  c.clear();
  c.resize(a.size() * b.size());

  size_t na = static_cast<size_t>(std::lround(std::sqrt(a.size())));
  assert(na * na == a.size());

  size_t nb = static_cast<size_t>(std::lround(std::sqrt(b.size())));
  assert(nb * nb == b.size());

  size_t nc = na * nb;
  assert(nc * nc == c.size());

  // [[0.3249 0.1083 0.1083 0.0361]
  //  [0.1083 0.0285 0.0361 0.0095]
  //  [0.1083 0.0361 0.0285 0.0095]
  //  [0.0361 0.0095 0.0095 0.0025]]

  for (size_t ar = 0; ar < na; ar++) {
    for (size_t ac = 0; ac < na; ac++) {
      for (size_t br = 0; br < nb; br++) {
        for (size_t bc = 0; bc < nb; bc++) {
          // std::cout << ar << " " << ac << " " << br << " " << bc <<
          // std::endl;
          c[(ar * nb + br) * nc + (ac * nb + bc)] =
              a[ar * na + ac] * b[br * nb + bc];
        }
      }
    }
  }
}

void get_probabilities(size_t n_vertices, std::vector<double>& probabilities) {
  std::vector<double> probabilities_base = {0.57, 0.19, 0.19, 0.05};
  probabilities = probabilities_base;
  {
    size_t scale = static_cast<size_t>(std::log2(n_vertices));
    for (size_t i = 1; i < scale; i++) {
      std::vector<double> prob_copy = probabilities;
      kronecker_prod(probabilities, prob_copy, probabilities_base);
    }
  }
}

double calculate_prob_error(std::vector<double>& empirical_prob) {
  size_t n = static_cast<size_t>(std::lround(std::sqrt(empirical_prob.size())));
  // Make sure we're a power of two
  assert(n > 0 && ((n & (n - 1)) == 0));

  std::vector<double> trusted_prob;
  get_probabilities(n, trusted_prob);

  std::cout << "Probabilities\n";
  print_matrix(trusted_prob);

  std::cout << "\nEstimated probabilities\n";
  print_matrix(empirical_prob);

  std::vector<double> abs_diff(empirical_prob.size(), 0);
  std::transform(empirical_prob.cbegin(), empirical_prob.cend(),
                 trusted_prob.cbegin(), abs_diff.begin(),
                 [&](double a, double b) { return std::abs(a - b); });

  std::cout << "\nDiff:\n";
  print_matrix(abs_diff);
  double l2_norm = std::sqrt(std::accumulate(
      abs_diff.begin(), abs_diff.end(), 0.0,
      [&](double prev, double next) { return prev + next * next; }));

  std::cout << "\nl^2-norm " << l2_norm << std::endl;

  return l2_norm;
}

#endif