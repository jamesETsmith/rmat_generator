#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest/doctest.h"
#include "rmat_generator.hpp"

TEST_CASE("initial") {
  rmat_generator my_gen(2);
  CHECK(my_gen.psum[0] == 0);
}