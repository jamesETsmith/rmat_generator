#include <algorithm>
#include <argparse/argparse.hpp>
#include <cassert>
#include <cmath>
#include <iomanip>  // std::setprecision
#include <iostream>
#include <mutex>
#include <thread>

#include "__rmat_utils.hpp"
#include "rmat_generator.hpp"

template <typename Generator>
void run_bench(size_t const n_vertices, size_t const n_threads,
               size_t const n_trials) {
  Generator base_gen(n_vertices);
  std::vector<std::thread> threads(n_threads);

  for (size_t i = 0; i < n_threads; i++) {
    threads[i] = std::thread([&]() {
      auto mygen = base_gen;
      mygen.set_seed(i);
      for (size_t i = 0; i < n_trials / n_threads; i++) {
        auto edge = mygen.next_edge();
        edge[0] += 1;  // Do something with edge so it's not unused
      }
    });
  }

  std::for_each(threads.begin(), threads.end(),
                std::mem_fn(&std::thread::join));
}

int main(int argc, char** argv) {
  argparse::ArgumentParser program(__FILE__);

  program.add_argument("n_vertices")
      .help("number of vertices in graph")
      .default_value(2)
      .scan<'i', int>();

  program.add_argument("n_trials")
      .help("number of trials to accumulate edges")
      .default_value(100000)
      .scan<'i', int>();

  program.add_argument("--n_threads")
      .help("number of threads used to generate edges")
      .default_value(12)
      .scan<'i', int>();

  program.add_argument("--engine")
      .help("PRNG to use. 0=std::mt19937, 1=std::ranlux48_base, 2=minstd_rand")
      .default_value(0)
      .scan<'i', int>();

  try {
    program.parse_args(argc, argv);
  } catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  size_t const n_vertices = program.get<int>("n_vertices");
  size_t const n_trials = program.get<int>("n_trials");
  size_t const n_threads = program.get<int>("--n_threads");

  enum class engine_opts {
    mt19937 = 0,
    ranlux48_base = 1,
    minstd_rand = 2,
    pcg = 3
  };
  engine_opts engine = static_cast<engine_opts>(program.get<int>("--engine"));

  switch (engine) {
    case engine_opts::mt19937: {
      std::cout << "Using std::mt19937 engine" << std::endl;
      run_bench<rmat_generator<std::mt19937>>(n_vertices, n_threads, n_trials);
      break;
    }
    case engine_opts::ranlux48_base: {
      std::cout << "Using std::ranlux48_base engine" << std::endl;
      run_bench<rmat_generator<std::ranlux48_base>>(n_vertices, n_threads,
                                                    n_trials);
      break;
    }
    case engine_opts::minstd_rand: {
      std::cout << "Using std::minstd_rand engine" << std::endl;
      run_bench<rmat_generator<std::minstd_rand>>(n_vertices, n_threads,
                                                  n_trials);
      break;
    }
    case engine_opts::pcg: {
      std::cout << "Using pcg32 engine" << std::endl;
      run_bench<rmat_generator<pcg32>>(n_vertices, n_threads, n_trials);
      break;
    }
  }
  // Error analysis
  //   calculate_prob_error(acc);

  return 0;
}