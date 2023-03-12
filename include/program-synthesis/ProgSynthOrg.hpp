#pragma once

#include <algorithm>
#include <tuple>

#include "emp/base/vector.hpp"
#include "TestResult.hpp"

// TODO - move genome and phenotype out of this file

namespace psynth {

template<typename PROGRAM_T>
struct ProgSynthGenome {
  using this_t = ProgSynthGenome<PROGRAM_T>;
  using program_t = PROGRAM_T;

  program_t program;

  ProgSynthGenome(const program_t& p) : program(p) {}
  ProgSynthGenome(const this_t&) = default;
  ProgSynthGenome(this_t&&) = default;

  bool operator==(const this_t& other) const {
    return program == other.program;
  }

  bool operator!=(const this_t& other) const {
    return !(*this == other);
  }

  bool operator<(const this_t& other) const {
    return program < other.program;
  }

  const program_t& GetProgram() const { return program; }
  program_t& GetProgram() { return program; }

};

// TODO - make this a class that better handles updates
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

};

template<typename PROGRAM_T>
class ProgSynthOrg {
public:
  using this_t = ProgSynthOrg<PROGRAM_T>;
  using program_t = PROGRAM_T;
  using phenotype_t = ProgSynthPhenotype;
  using genome_t = ProgSynthGenome<PROGRAM_T>;

protected:
  phenotype_t phenotype;
  genome_t genome;

  size_t pop_id = 0;

public:

  ProgSynthOrg(const genome_t& g) :
    phenotype(),
    genome(g)
  { ; }

  ProgSynthOrg(const ProgSynthOrg&) = default;
  ProgSynthOrg(ProgSynthOrg&&) = default;

  genome_t& GetGenome() { return genome; }
  const genome_t& GetGenome() const { return genome; }

  // phenotype_t& GetPhenotype() { return phenotype; }
  const phenotype_t& GetPhenotype() const { return phenotype; }

  size_t GetPopID() const { return pop_id; }
  void SetPopID(size_t id) { pop_id = id; }

  void ResetPhenotype(size_t num_tests=1) {
    phenotype.Reset(num_tests);
  }

  void UpdatePhenotype(size_t test_id, const TestResult& test_result) {
    phenotype.Update(test_id, test_result);
  }

};

}