#pragma once

#include <algorithm>
#include "BaseDiagnostic.hpp"

namespace diag {

struct MultiPathExplorationDiagnostic : public BaseDiagnostic {
  using phenotype_t = emp::vector<double>;
  using genotype_t = emp::vector<double>;

  phenotype_t Translate(const genotype_t& genome) const {
    phenotype_t phenotype(genome.size(), 0);
    Translate(genome, phenotype);
    return phenotype;
  }

  void Translate(const genotype_t& genome, phenotype_t& phenotype) const {
    // phenotype should be same size as genome
    phenotype.resize(genome.size(), 0);
    // fill phenotype with baseline credit value
    const double fill_value = GetMaxErrorCredit();
    std::fill(
      phenotype.begin(),
      phenotype.end(),
      fill_value
    );
    // find max value position
    const auto max_gene_it = std::max_element(genome.begin(), genome.end());
    const size_t max_gene = std::distance(
      genome.begin(),
      max_gene_it
    );

    // find where the order breaks
    const size_t sort = std::distance(
      genome.begin(),
      std::is_sorted_until(
        max_gene_it,
        genome.end(),
        std::greater<>()
      )
    );
    // Middle of optimal value till order broken
    for (size_t i = max_gene; i < sort; ++i) {
      phenotype[i] = genome[i];
    }
  }

};

}