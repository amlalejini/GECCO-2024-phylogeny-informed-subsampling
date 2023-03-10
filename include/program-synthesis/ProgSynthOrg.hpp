#pragma once

#include <algorithm>
#include <tuple>

#include "emp/base/vector.hpp"

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
struct ProgSynthPhenotype {
  using this_t = ProgSynthPhenotype;

  // TODO - replace three vectors with one vector of TestResults?
  emp::vector<double> test_scores;   ///< Scores on tests.
  emp::vector<bool> test_passes;     ///< Pass/fail on tests
  emp::vector<bool> test_evaluated;  ///< Whether a test has been evaluated.

  double aggregate_score = 0.0;
  size_t num_passes = 0;
  bool is_solution = false;

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

  bool operator==(const this_t& o) const {
    return std::tie(
      aggregate_score,
      test_scores,
      test_evaluated,
      num_passes,
      is_solution
    ) == std::tie(
      o.aggregate_score,
      o.test_scores,
      o.test_evaluated,
      o.num_passes,
      o.is_solution
    );
  }

  bool operator!=(const this_t& o) const {
    return !(*this == o);
  }

  bool operator<(const this_t& o) const {
    return std::tie(
      aggregate_score,
      test_scores,
      test_evaluated,
      num_passes,
      is_solution
    ) < std::tie(
      o.aggregate_score,
      o.test_scores,
      o.test_evaluated,
      o.num_passes,
      o.is_solution
    );
  }

  double GetAggregateScore() const {
    return aggregate_score;
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

  bool evaluated = false;

  size_t pop_id = 0;

public:

  ProgSynthOrg(const genome_t& g) :
    phenotype(),
    genome(g)
  {
    evaluated = false;
  }

  ProgSynthOrg(const ProgSynthOrg&) = default;
  ProgSynthOrg(ProgSynthOrg&&) = default;

  genome_t& GetGenome() { return genome; }
  const genome_t& GetGenome() const { return genome; }

  phenotype_t& GetPhenotype() { return phenotype; }
  const phenotype_t& GetPhenotype() const { return phenotype; }

  size_t GetPopID() const { return pop_id; }
  void SetPopID(size_t id) { pop_id = id; }

  bool IsEvaluated() const { return evaluated; }

  // TODO - process phenotype
};

}