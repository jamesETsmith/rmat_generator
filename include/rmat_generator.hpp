#ifndef __RMAT_GENERTAOR_HPP
#define __RMAT_GENERTAOR_HPP

// #include <algorithm>
#include <array>
// #include <numeric>
#include <cmath>
#include <random>

class rmat_generator {
  //  private:
 public:
  size_t const n_vertices;
  double const a, b, c, d;
  std::array<double, 4> const psum;
  std::mt19937 gen;
  std::uniform_real_distribution<> dis{0.0, 1.0};

  //  public:
  rmat_generator(size_t n_vertices, double a = 0.57, double b = 0.19,
                 double c = 0.19, double d = 0.05,
                 size_t seed = std::random_device{}())
      : n_vertices(n_vertices),
        a(a),
        b(b),
        c(c),
        d(d),
        // gen(std::mt19937(seed)),
        psum(std::array<double, 4>{0, a, a + b, a + b + c}) {
    gen.seed(seed);
  };

  std::array<size_t, 2> next_edge() {
    size_t rl = 0, ru = n_vertices;
    size_t cl = 0, cu = n_vertices;

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