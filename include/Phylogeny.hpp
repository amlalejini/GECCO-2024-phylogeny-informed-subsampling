#pragma once

#include <map>
#include <set>
#include <unordered_set>

#include "emp/Evolve/Systematics.hpp"
#include "emp/base/vector.hpp"

// TODO - Data structure for tracking phenotypes
// TODO - Subclass of systematics that makes it easy to search for relative

// NOTE - might need to make this more generic!
struct phenotype_info {

  using phen_t = emp::vector<double>;
  using has_phen_t = std::true_type;
  using has_mutations_t = std::false_type;
  using has_fitness_t = std::true_type;

  // Phenotype is set of test case scores
  // Need to set an eval vector

  emp::DataNode<
    double,
    emp::data::Current,
    emp::data::Range
  > fitness;            ///< Tracks fitness range for this taxon.

  phen_t phenotype;     ///< Tracks phenotype associated with this taxon.
  emp::vector<bool> traits_evaluated; ///< Tracks whether a particular trait has been evaluated

  const phen_t& GetPhenotype() const {
    return phenotype;
  }

  double GetFitness() const {
    return fitness.GetMean();
  }

  const emp::vector<bool>& GetTraitsEvaluated() const {
    return traits_evaluated;
  }

  bool TraitEvaluated(size_t i) const {
    emp_assert(i < traits_evaluated.size());
    return traits_evaluated[i];
  }

  void RecordFitness(double fit) {
    fitness.Add(fit);
  }

  // TODO - test if this makes a proper copy
  void RecordPhenotype(
    const phen_t& phen,
    const emp::vector<bool>& eval
  ) {
    phenotype = phen;
    traits_evaluated = eval;
  }

};