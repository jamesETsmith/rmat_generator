#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "__rmat_utils.hpp"
#include "doctest/doctest.h"
#include "rmat_generator.hpp"

TEST_CASE("n_vertices=2") {
  size_t const seed = 2;
  size_t const n_vertices = 2;
  size_t const n_trials = 2000000;
  rmat_generator mygen(n_vertices, seed);

  std::vector<double> acc(n_vertices * n_vertices);
  for (size_t i = 0; i < n_trials; i++) {
    auto edge = mygen.next_edge();
    acc[edge[0] * n_vertices + edge[1]] += 1. / static_cast<double>(n_trials);
  }

  double l2_norm = calculate_prob_error(acc);
  CHECK(l2_norm < 1e-3);
}

TEST_CASE("n_vertices=4") {
  size_t const n_vertices = 4;
  size_t const seed = n_vertices;
  size_t const n_trials = 10000000;
  rmat_generator mygen(n_vertices, seed);

  std::vector<double> acc(n_vertices * n_vertices);
  for (size_t i = 0; i < n_trials; i++) {
    auto edge = mygen.next_edge();
    acc[edge[0] * n_vertices + edge[1]] += 1. / static_cast<double>(n_trials);
  }

  double l2_norm = calculate_prob_error(acc);
  CHECK(l2_norm < 1e-3);
}

TEST_CASE("n_vertices=8") {
  size_t const n_vertices = 8;
  size_t const seed = n_vertices;
  size_t const n_trials = 10000000;
  rmat_generator mygen(n_vertices, seed);

  std::vector<double> acc(n_vertices * n_vertices);
  for (size_t i = 0; i < n_trials; i++) {
    auto edge = mygen.next_edge();
    acc[edge[0] * n_vertices + edge[1]] += 1. / static_cast<double>(n_trials);
  }

  double l2_norm = calculate_prob_error(acc);
  CHECK(l2_norm < 1e-3);
}