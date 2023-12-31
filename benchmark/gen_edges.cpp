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

  try {
    program.parse_args(argc, argv);
  } catch (const std::exception& err) {
    std::cerr << err.what() << std::endl;
    std::cerr << program;
    return 1;
  }

  size_t n_vertices = static_cast<size_t>(program.get<int>("n_vertices"));
  size_t n_trials = static_cast<size_t>(program.get<int>("n_trials"));

  std::vector<double> acc(n_vertices * n_vertices, 0);
  std::mutex m_acc;

  size_t const N = static_cast<size_t>(program.get<int>("--n_threads"));
  std::vector<std::thread> threads(N);

  for (size_t i = 0; i < N; i++) {
    threads[i] = std::thread([&]() {
      std::vector<double> local_acc(acc.size(), 0);
      rmat_generator mygen(n_vertices, i);

      for (size_t i = 0; i < n_trials / N; i++) {
        auto edge = mygen.next_edge();
        local_acc[edge[0] * n_vertices + edge[1]] +=
            1. / static_cast<double>(n_trials);
      }

      // Updated the shared accumulator
      std::lock_guard<std::mutex> lock(m_acc);
      std::transform(local_acc.cbegin(), local_acc.cend(), acc.cbegin(),
                     acc.begin(), [](double a, double b) { return a + b; });
    });
  }

  std::for_each(threads.begin(), threads.end(),
                std::mem_fn(&std::thread::join));

  // Error analysis
  calculate_prob_error(acc);

  return 0;
}