/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/
#pragma once

#include <algorithm>
#include "BaseDiagnostic.hpp"

namespace diag {

struct ContradictoryObjectivesDiagnostic : public BaseDiagnostic {
  using phenotype_t = emp::vector<double>;
  using genotype_t = emp::vector<double>;

  ContradictoryObjectivesDiagnostic(double targ, double err_cred)
    : BaseDiagnostic(targ, err_cred)
  { ; }
  ContradictoryObjectivesDiagnostic() = default;

  phenotype_t Translate(const genotype_t& genome) const {
    phenotype_t phenotype(genome.size(), 0);
    Translate(genome, phenotype);
    return phenotype;
  }

  void Translate(const genotype_t& genome, phenotype_t& phenotype) const {
    // phenotype should be same size as genome
    phenotype.resize(genome.size(), 0);
    const double fill_value = GetMaxErrorCredit();
    std::fill(
      phenotype.begin(),
      phenotype.end(),
      fill_value
    );
    // find max value position
    const size_t max_gene = std::distance(
      genome.begin(),
      std::max_element(genome.begin(), genome.end())
    );
    // Set max value position to gene value
    phenotype[max_gene] = genome[max_gene];
  }

};

}