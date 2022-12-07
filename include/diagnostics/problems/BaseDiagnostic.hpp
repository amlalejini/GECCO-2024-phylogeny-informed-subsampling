/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/
#pragma once

#include <algorithm>
#include "emp/base/vector.hpp"

namespace diag {

  struct BaseDiagnostic {
    using phenotype_t = emp::vector<double>;
    using genotype_t = emp::vector<double>;

    double target = 0.0;           ///< Target trait value.
    double max_error_credit = 0.0; ///< Maximum credit allowed for error.

    BaseDiagnostic(double targ, double err_cred)
      : target(targ),
        max_error_credit(err_cred)
    { ; }

    BaseDiagnostic() = default;
    ~BaseDiagnostic() = default;

    double GetTarget() const { return target; }
    double GetMaxErrorCredit() const { return max_error_credit; }

    void SetTarget(double t) {
      target = t;
    }

    void SetMaxErrorCredit(double c) {
      max_error_credit = c;
    }

    virtual phenotype_t Translate(const genotype_t& genome) const = 0;
    virtual void Translate(const genotype_t& genome, phenotype_t& phenotype) const = 0;

  };

}