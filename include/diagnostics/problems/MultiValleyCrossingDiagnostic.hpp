#pragma once

#include <algorithm>
#include "BaseDiagnostic.hpp"

namespace diag {

struct MultiValleyCrossingDiagnostic : public BaseDiagnostic {
  using phenotype_t = emp::vector<double>;
  using genotype_t = emp::vector<double>;

  // Defaults used by Hernandez et al
  double dips_start = 8.0;
  double dips_end = 99.9;
  emp::vector<double> peaks = {
    -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0, -1.0,  8.0,  9.0,
    9.0, 11.0, 11.0, 11.0, 14.0, 14.0, 14.0, 14.0, 18.0, 18.0,
    18.0, 18.0, 18.0, 23.0, 23.0, 23.0, 23.0, 23.0, 23.0, 29.0,
    29.0, 29.0, 29.0, 29.0, 29.0, 29.0, 36.0, 36.0, 36.0, 36.0,
    36.0, 36.0, 36.0, 36.0, 44.0, 44.0, 44.0, 44.0, 44.0, 44.0,
    44.0, 44.0, 44.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0, 53.0,
    53.0, 53.0, 53.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0,
    63.0, 63.0, 63.0, 63.0, 74.0, 74.0, 74.0, 74.0, 74.0, 74.0,
    74.0, 74.0, 74.0, 74.0, 74.0, 74.0, 86.0, 86.0, 86.0, 86.0,
    86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 86.0, 99.0
  };

  MultiValleyCrossingDiagnostic(
    double targ,
    double err_cred,
    double _dips_start,
    double _dips_end,
    const emp::vector<double>& _peaks
  ) :
    BaseDiagnostic(targ, err_cred),
    dips_start(_dips_start),
    dips_end(_dips_end),
    peaks(_peaks)
  { ; }

  MultiValleyCrossingDiagnostic(
    double targ,
    double err_cred
  ) :
    BaseDiagnostic(targ, err_cred)
  { ; }

  MultiValleyCrossingDiagnostic() = default;

  phenotype_t Translate(const genotype_t& genome) const {
    emp_assert(genome.size() == peaks.size(), "Genome size must match configured peaks size.");
    phenotype_t phenotype(genome.size(), 0);
    Translate(genome, phenotype);
    return phenotype;
  }

  void Translate(const genotype_t& genome, phenotype_t& phenotype) const {
    emp_assert(genome.size() == peaks.size(), "Genome size must match configured peaks size.");
    // phenotype should be same size as genome
    phenotype.resize(genome.size(), 0);
    // fill phenotype with baseline credit value
    const double fill_value = GetMaxErrorCredit();
    std::fill(
      phenotype.begin(),
      phenotype.end(),
      fill_value
    );

    for (size_t i = 0; i < phenotype.size(); ++i) {
      phenotype[i] = (phenotype[i] <= dips_start || phenotype[i] >= dips_end)
        ? phenotype[i]
        : 2.0 * peaks[static_cast<size_t>(phenotype[i])] - phenotype[i];
    }
  }

};

}