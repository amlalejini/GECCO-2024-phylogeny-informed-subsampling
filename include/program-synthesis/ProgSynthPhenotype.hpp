#pragma once

#include <utility>
#include <algorithm>

#include "emp/base/vector.hpp"

#include "TestResult.hpp"

namespace psynth {

class ProgSynthPhenotype {
public:
  using this_t = ProgSynthPhenotype;
protected:

  emp::vector<double> test_scores;   ///< Scores on tests.
  emp::vector<bool> test_passes;     ///< Pass/fail on tests
  emp::vector<bool> test_evaluated;  ///< Whether a test has been evaluated.

  double aggregate_score = 0.0;
  size_t num_passes = 0;
  bool is_solution = false;

public:
  void Reset(size_t num_tests=1) {
    aggregate_score = 0.0;
    num_passes = 0;
    is_solution = false;

    test_scores.resize(num_tests);
    std::fill(
      test_scores.begin(),
      test_scores.end(),
      0
    );

    test_evaluated.resize(num_tests);
    std::fill(
      test_evaluated.begin(),
      test_evaluated.end(),
      false
    );

    test_passes.resize(num_tests);
    std::fill(
      test_passes.begin(),
      test_passes.end(),
      false
    );
  }

  void Update(size_t test_id, const TestResult& test_result, bool evaluated=true) {
    test_scores[test_id] = test_result.score;
    test_passes[test_id] = test_result.is_correct;
    test_evaluated[test_id] = evaluated;
    num_passes += (size_t)test_result.is_correct;
    aggregate_score += test_result.score;
  }

  bool operator==(const this_t& o) const {
    // Only use scores and is solution to test equality
    return std::tie(
      test_scores,
      is_solution
    ) == std::tie(
      o.test_scores,
      o.is_solution
    );
  }

  bool operator!=(const this_t& o) const {
    return !(*this == o);
  }

  bool operator<(const this_t& o) const {
    // Only use scores and is solution for comparison
    return std::tie(
      test_scores,
      is_solution
    ) < std::tie(
      o.test_scores,
      o.is_solution
    );
  }

  double GetAggregateScore() const {
    return aggregate_score;
  }

  bool PassedTest(size_t test_id) const {
    return test_passes[test_id];
  }

  double GetTestScore(size_t test_id) const {
    return test_scores[test_id];
  }

  bool TestEvaluated(size_t test_id) const {
    return test_evaluated[test_id];
  }

  const emp::vector<double>& GetTestScores() const {
    return test_scores;
  }

};

}