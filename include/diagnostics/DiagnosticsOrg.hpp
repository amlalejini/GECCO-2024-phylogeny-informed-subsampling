/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/
#pragma once

#include <algorithm>
#include <numeric>
#include <functional>
#include <unordered_map>
#include "emp/base/vector.hpp"

namespace diag {

class DiagnosticsOrg {
public:
  using genome_t = emp::vector<double>;
  using phenotype_t = emp::vector<double>;

protected:
  genome_t genome;
  // can only modify phenotype via SetPhenotype and TranslateGenome functions
  phenotype_t phenotype;
  bool evaluated = false;

  emp::vector<bool> optimal_traits;
  bool optimal_traits_calculated=false;
  size_t streak_size = 0;

  double aggregate_score = 0.0;
  size_t max_trait_id = 0;

  size_t pop_id = 0;

  std::unordered_map<std::string, size_t> muts_from_parent;

  // bool is_clone = false;

  void AggregatePhenotype() {
    aggregate_score = std::accumulate(
      phenotype.begin(),
      phenotype.end(),
      0.0
    );
  }

  size_t FindMaxTraitID() {
    return std::distance(
      phenotype.begin(),
      std::max_element(phenotype.begin(), phenotype.end())
    );
  }

  size_t FindStreakSize() {
    // get longest streak
    size_t count = 0;
    size_t max_cnt = 0;
    for (auto s : phenotype) {
      if (s > 0.0) {
        count++;
      } else {
        if (count > max_cnt) {
          max_cnt = count;
        }
        count = 0;
      }
    }
    if (count > max_cnt) {
      max_cnt = count;
    }
    return max_cnt;
  }

public:

  DiagnosticsOrg(size_t num_genes, double init_gene_val=0.0) {
    genome.resize(num_genes, init_gene_val);
    phenotype.resize(num_genes, 0);
    optimal_traits.resize(num_genes, false);
    evaluated = false;
  }

  DiagnosticsOrg(const genome_t& g)
    : genome(g)
  {
    phenotype.resize(genome.size(), 0);
    optimal_traits.resize(genome.size(), false);
    evaluated = false;
  }

  DiagnosticsOrg(const DiagnosticsOrg&) = default;
  DiagnosticsOrg(DiagnosticsOrg&&) = default;
  ~DiagnosticsOrg() { ; }
  DiagnosticsOrg& operator=(const DiagnosticsOrg&) = default;
  DiagnosticsOrg& operator=(DiagnosticsOrg&&) = default;

  const genome_t& GetGenome() const {
    return genome;
  }

  genome_t & GetGenome() {
    emp_assert(genome.size() > 0);
    return genome;
  }

  const phenotype_t& GetPhenotype() const {
    return phenotype;
  }

  size_t GetPopID() const { return pop_id; }
  void SetPopID(size_t id) { pop_id = id; }

  bool IsEvaluated() const { return evaluated; }

  const emp::vector<bool>& GetOptimalTraits() const {
    emp_assert(evaluated);
    return optimal_traits;
  }

  bool IsGeneOptimal(size_t id) const {
    emp_assert(evaluated);
    emp_assert(id < optimal_traits.size());
    return optimal_traits[id];
  }

  double GetAggregateScore() const {
    emp_assert(evaluated);
    return aggregate_score;
  }

  size_t GetStreakSize() const {
    emp_assert(evaluated);
    return streak_size;
  }

  size_t GetMaxTraitID() const {
    emp_assert(evaluated);
    return max_trait_id;
  }

  // bool IsClone() const {
  //   return is_clone;
  // }

  const std::unordered_map<std::string, size_t>& GetMutsFromParent() const {
    return muts_from_parent;
  }

  size_t GetMutsFromParent(const std::string& type) const {
    emp_assert(emp::Has(muts_from_parent, type));
    return muts_from_parent.at(type);
  }

  void SetMutsFromParent(const std::string& type, size_t count) {
    muts_from_parent[type] = count;
  }

  void SetPhenotype(const phenotype_t& phen) {
    emp_assert(phen.size() == genome.size());
    // Update phenotype
    phenotype = phen;
    // aggregate score
    AggregatePhenotype();
    // find id of max gene
    max_trait_id = FindMaxTraitID();
    // mark evaluated to be true
    evaluated = true;
  }

  void TranslateGenome(
    const std::function<void(const genome_t&, phenotype_t&)>& translate
  ) {
    translate(genome, phenotype);
    // aggregate score
    AggregatePhenotype();
    // find id of max gene
    max_trait_id = FindMaxTraitID();
    // mark evaluated to be true
    evaluated = true;
  }

  void CalcOptimalTraits(double target, double accuracy) {
    emp_assert(evaluated);
    optimal_traits.resize(phenotype.size(), false);
    // Identify traits within accuracy threshold of target
    emp_assert(target > 0, target);
    emp_assert(accuracy > 0 && accuracy <= 1.0, accuracy, "Accuracy must be between 0 and 1");
    const double threshold = accuracy * target;
    for (size_t i = 0; i < phenotype.size(); ++i) {
      optimal_traits[i] = phenotype[i] >= threshold;
    }
    streak_size = FindStreakSize();
    optimal_traits_calculated = true;
  }

};

}

