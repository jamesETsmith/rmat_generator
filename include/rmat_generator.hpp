#ifndef __RMAT_GENERTAOR_HPP
#define __RMAT_GENERTAOR_HPP

#include <array>
#include <cassert>
#include <cmath>
#include <limits>
#include <random>

// For compatability
#ifdef __EMU_CC__
#undef noinline
#define noinline __noinline__
#endif

#include "pcg_random.hpp"

namespace detail {

/**
 * @brief Return the number of bits set to 1 in n
 *
 * @param n
 * @return size_t
 */
constexpr inline size_t popcount(size_t n) {
  size_t count = 0;
  while (n) {
    count += n & 1;
    n >>= 1;
  }
  return count;
}

/**
 * @brief Calculate log2 for a size_t (either 32- or 64-bit). Doesn't handle the
 * case when N is 0 or warn the users if the log2(N) is non-integral.
 *
 * @param N
 * @return size_t
 */
constexpr inline size_t fast_log2_unsafe(size_t N) {
  N |= (N >> 1);
  N |= (N >> 2);
  N |= (N >> 4);
  N |= (N >> 8);
  N |= (N >> 16);

  if constexpr (sizeof(size_t) == 8) {
    N |= (N >> 32);
  }

  return popcount(N) - 1;
}

/**
 * @brief Return 2^N
 *
 * @param N
 * @return size_t
 */
constexpr inline size_t pow_2(size_t N) { return 1ULL << N; }

/**
 * @brief Ripped unceremoniously from LLVM's libc++ and modified for several
 * reasons:
 *
 * 1) Add constexpr where possible to make sure we're calculating the
 * minimum amount at runtime.
 *
 * 2) GCC's implementation of uniform_real_dist has several calls to std::log
 * with are NOT constexpr.
 *
 * See
 * https://github.com/llvm/llvm-project/blob/main/libcxx/include/__random/uniform_real_distribution.h
 *
 * _URNG must look and behave like the random number generators in STL. (NO
 * CHECKS FOR THIS NOW)
 *
 * @tparam _RealType
 * @tparam _URNG
 * @param __g
 * @return _RealType
 */
template <class _RealType, class _URNG>
inline _RealType stripped_uniform_real_dist(_URNG& __g) {
  constexpr size_t __dt = std::numeric_limits<_RealType>::digits;
  // Compute the number of bits needed to represent one random value from the
  // generator
  constexpr size_t __log_r =
      fast_log2_unsafe(_URNG::max() - _URNG::min() + uint64_t(1));
  // __log2<uint64_t, _URNG::max() - _URNG::min() + uint64_t(1)>::value;

  // Determine the number of random values needed to generate __dt bits
  constexpr const size_t __k =
      __dt / __log_r + (__dt % __log_r != 0) + (__dt == 0);

  // Compute the range of the random number generator
  constexpr const _RealType __rp =
      static_cast<_RealType>(_URNG::max() - _URNG::min()) + _RealType(1);

  _RealType __base = __rp;
  _RealType __sp = __g() - _URNG::min();
  for (size_t __i = 1; __i < __k; ++__i, __base *= __rp)
    __sp += (__g() - _URNG::min()) * __base;
  return __sp / __base;
}

};  // namespace detail

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
  size_t const scale;
  double const a, b, c, d;

  std::array<double, 4> const psum;
  Generator gen;

 public:
  rmat_generator(size_t n_vertices, size_t seed = std::random_device{}(),
                 double a = 0.57, double b = 0.19, double c = 0.19,
                 double d = 0.05)
      : n_vertices(n_vertices),
        scale(detail::fast_log2_unsafe(n_vertices)),
        a(a),
        b(b),
        c(c),
        d(d),
        psum(std::array<double, 4>{0, a, a + b, a + b + c}) {
    // Input checking
    // n_vertices must be a power of two
    assert(n_vertices == detail::pow_2(scale));

    // Probabilities must be 1
    // watch out for the epsilon
    assert(std::abs(a + b + c + d - 1.0) <
           std::numeric_limits<double>::epsilon());

    // Set seed
    gen.seed(seed);
  }

  void set_seed(size_t seed) { gen.seed(seed); }

  std::array<size_t, 2> next_edge() {
    size_t dst = 0, src = 0;
    for (size_t bit = ((size_t)1) << (scale - 1); bit >= 1; bit >>= 1) {
      double r = detail::stripped_uniform_real_dist<double>(gen);
      if (r > a) {      /* outside quadrant 1 */
        if (r <= a + b) /* in quadrant 2 */
          dst |= bit;
        else if (r <= a + b + c) /* in quadrant 3 */
          src |= bit;
        else { /* in quadrant 4 */
          dst |= bit;
          src |= bit;
        }
      }
    }
    return {src, dst};
  }

  std::array<size_t, 2> next_edge_noisy() {
    size_t dst = 0, src = 0;
    double A = a, B = b, C = c, D = d;
    for (size_t bit = ((size_t)1) << (scale - 1); bit >= 1; bit >>= 1) {
      double r = detail::stripped_uniform_real_dist<double>(gen);
      if (r > A) {      /* outside quadrant 1 */
        if (r <= A + B) /* in quadrant 2 */
          dst |= bit;
        else if (r <= A + B + C) /* in quadrant 3 */
          src |= bit;
        else { /* in quadrant 4 */
          dst |= bit;
          src |= bit;
        }
      }

      /*
        Assuming R is in [0, 1), 0.9 + 0.2 * R is in (0.90, 1.10).
        So the new probabilities are the old +/- 10%
      */
      A *= (9. + 2 * detail::stripped_uniform_real_dist<double>(gen)) / 10;
      B *= (9. + 2 * detail::stripped_uniform_real_dist<double>(gen)) / 10;
      C *= (9. + 2 * detail::stripped_uniform_real_dist<double>(gen)) / 10;
      D *= (9. + 2 * detail::stripped_uniform_real_dist<double>(gen)) / 10;
      /* Used 5 random numbers. */

      {
        const double norm = 1.0 / (A + B + C + D);
        A *= norm;
        B *= norm;
        C *= norm;
      }
      /* So long as +/- are monotonic, ensure a+b+c+d <= 1.0 */
      D = 1.0 - (A + B + C);
    }
    return {src, dst};
  }
};

#endif