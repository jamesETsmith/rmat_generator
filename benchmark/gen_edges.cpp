#include <algorithm>
#include <argparse/argparse.hpp>
#include <cassert>
#include <cmath>
#include <iomanip>  // std::setprecision
#include <iostream>
#include <mutex>
#include <thread>

#include "rmat_generator.hpp"

void print_matrix(std::vector<double>& mat) {
  size_t n = std::lround(std::sqrt(mat.size()));
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

  size_t na = std::lround(std::sqrt(a.size()));
  assert(na * na == a.size());

  size_t nb = std::lround(std::sqrt(b.size()));
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

  size_t n_vertices = program.get<int>("n_vertices");
  size_t n_trials = program.get<int>("n_trials");

  std::vector<double> acc(n_vertices * n_vertices, 0);
  std::mutex m_acc;

  size_t const N = program.get<int>("--n_threads");
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

  std::vector<double> probabilities_base = {0.57, 0.19, 0.19, 0.05};
  std::vector<double> probabilities = probabilities_base;
  {
    size_t n_iter = static_cast<size_t>(std::log2(n_vertices));
    for (size_t i = 1; i < n_iter; i++) {
      std::vector<double> prob_copy = probabilities;
      kronecker_prod(probabilities, prob_copy, probabilities_base);
    }
  }

  std::cout << "Probabilities\n";
  print_matrix(probabilities);

  std::cout << "\nEstimated probabilities\n";
  print_matrix(acc);

  std::vector<double> abs_diff(acc.size(), 0);
  std::transform(acc.cbegin(), acc.cend(), probabilities.cbegin(),
                 abs_diff.begin(),
                 [&](double a, double b) { return std::abs(a - b); });

  std::cout << "\nDiff:\n";
  print_matrix(abs_diff);
  double l2_norm = std::sqrt(std::accumulate(
      abs_diff.begin(), abs_diff.end(), 0.0,
      [&](double prev, double next) { return prev + next * next; }));

  std::cout << "\nl^2-norm " << l2_norm << std::endl;
  return 0;
}