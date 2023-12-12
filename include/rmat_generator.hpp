#ifndef __RMAT_GENERTAOR_HPP
#define __RMAT_GENERTAOR_HPP

// #include <algorithm>
#include <array>
// #include <numeric>
#include <cassert>
#include <cmath>
#include <limits>
#include <random>

#include "pcg_random.hpp"

/*
FROM: https://en.cppreference.com/w/cpp/numeric/random

The choice of which engine to use involves a number of trade-offs: the linear
congruential engine is moderately fast and has a very small storage requirement
for state. The lagged Fibonacci generators are very fast even on processors
without advanced arithmetic instruction sets, at the expense of greater state
storage and sometimes less desirable spectral characteristics. The Mersenne
twister is slower and has greater state storage requirements but with the right
parameters has the longest non-repeating sequence with the most desirable
spectral characteristics (for a given definition of desirable).
*/

template <typename Generator = pcg32>
class rmat_generator {
 private:
  size_t const n_vertices;
  double const a, b, c, d;
  std::array<double, 4> const psum;
  Generator gen;
  std::uniform_real_distribution<> dis{0.0, 1.0};

 public:
  rmat_generator(size_t n_vertices, size_t seed = std::random_device{}(),
                 double a = 0.57, double b = 0.19, double c = 0.19,
                 double d = 0.05)
      : n_vertices(n_vertices),
        a(a),
        b(b),
        c(c),
        d(d),
        psum(std::array<double, 4>{0, a, a + b, a + b + c}) {
    // Input checking
    // n_vertices must be a power of two
    size_t scale = static_cast<size_t>(std::log2(n_vertices));
    // try to use bit manipulation for this
    assert(static_cast<double>(n_vertices) == std::pow(2.0, scale));

    // Probabilities must be 1
    // watch out for the epsilon
    assert(abs(a + b + c + d - 1.0) < std::numeric_limits<double>::epsilon());

    // Set seed
    gen.seed(seed);
  };

  void set_seed(size_t seed) { gen.seed(seed); }

  inline std::array<size_t, 2> next_edge() {
    size_t rl = 0, ru = n_vertices;
    size_t cl = 0, cu = n_vertices;

    // switch to bit manipulation to calculate this
    size_t n_iter = static_cast<size_t>(std::log2(n_vertices));

    for (size_t i = 0; i < n_iter; i++) {
      double rn = dis(gen);
      if (rn < psum[1]) {
        ru = rl + (ru - rl) / 2;
        cu = cl + (cu - cl) / 2;
      } else if (rn < psum[2]) {
        ru = rl + (ru - rl) / 2;
        cl = cl + (cu - cl) / 2;
      } else if (rn < psum[3]) {
        rl = rl + (ru - rl) / 2;
        cu = cl + (cu - cl) / 2;
      } else {
        rl = rl + (ru - rl) / 2;
        cl = cl + (cu - cl) / 2;
      }
    }

    return {rl, cl};
  };
};

#endif