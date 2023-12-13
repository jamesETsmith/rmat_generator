#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "rmat_generator.hpp"

TEST_CASE("popcount") {
  size_t x1 = 0b01;

  CHECK(x1 == 1);
  CHECK(detail::popcount(x1) == 1);

  size_t x2 = 0b1111000;
  CHECK(detail::popcount(x2) == 4);

  size_t x3 = 0b00111000;
  CHECK(detail::popcount(x3) == 3);
}

TEST_CASE("fast_log2") {
  size_t x1 = 0b10;
  CHECK(detail::fast_log2_unsafe(x1) == 1);

  size_t x2 = 0b100;
  CHECK(detail::fast_log2_unsafe(x2) == 2);

  size_t x3 = 0b100000000000;
  CHECK(detail::fast_log2_unsafe(x3) == 11);

  size_t x4 = 0b0001;
  CHECK(detail::fast_log2_unsafe(x4) == 0);
}

TEST_CASE("pow_2") {
  size_t x1 = 2;
  CHECK(detail::pow_2(x1) == 4);

  size_t x2 = 1;
  CHECK(detail::pow_2(x2) == 2);

  size_t x3 = 16;
  CHECK(detail::pow_2(x3) == 65536);
}