/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/

#ifndef DIAGNOSTICS_PROBLEM_HPP
#define DIAGNOSTICS_PROBLEM_HPP

#include <algorithm>
#include "emp/base/vector.hpp"

namespace diag {

class Diagnostic {
public:
  using phenotype_t = emp::vector<double>;
  using score_t = emp::vector<double>;
  using genome_t = emp::vector<double>;

protected:
  // holds the target objective value
  double target = 0.0;
  // holds maximum credit allowed for error
  double max_cred;
  // max credit set?
  bool cred_set = false;

public:

  Diagnostic(double targ, double cred)
    : target(targ),
      max_cred(cred),
      cred_set(true)
  { ; }

  Diagnostic() = default;

  // getters
  double GetTarget() const { return target; }
  double GetMaxCredit() const { return max_cred; }

  // setters
  void SetTarget(double t) {
    target = t;
  }

  void SetMaxCredit(double c) {
    cred_set = true;
    max_cred = c;
  }


  /**
   * Exploitation function:
   *
   * All genes in genome evaluated treated as independent optimization tasks.
   * In other words, there is no interactions between genes.
   * The phenotype vector is just the genome, with the caveat that genes cannot go over associated target value.
   *
   * @param g genome from organism being evaluated.
   *
   * @return phenotype vector that is calculated from 'g'.
   */
  phenotype_t ExploitationRate(const genome_t& g);

  /**
   * Ordered Exploitation function:
   *
   * All genes in genome expected to follow descending order for evaluation.
   * While in decending order, phenotype[i] = genome[i].
   * If descending order is broken, max_credit is assigned to each following trait .
   * Note that this problem explicitly forces optimization from start to end.
   *
   * @param g genome from organism being evaluated.
   *
   * @return phenotype vector that is calculated from 'g'.
   */
  phenotype_t OrderedExploitation(const genome_t& g);

  /**
   * Contradictory Objectives function:
   *
   * Solutions are pressured to optimize one gene, while all other genes are ignored and set to max_credit.
   * To calculate phenotype vector, first we find the maximum gene value at position i.
   * After, phenotype[i] = genome[i], and every every other trait is set to max_credit.
   *
   * @param g genome from organism being evaluated.
   *
   * @return phenotype vector that is calculated from 'g'.
   */
  phenotype_t ContradictoryObjectives(const genome_t& g);

  /**
   * Multi-path Exploration function:
   *
   * Solutions are must be able to simultaneously explore multiple pathways, with the goal of pursuing the one that leads to the optimum.
   * To calculate phenotype vector, first we find the maximum gene value at position i.
   * From that gene forward, while genes follow decending order, phenotype[i] = genome[i].
   * This process stops until the decending order is broken, or the end of the genome is reached.
   * Note, the optimum is reached when i = 0, the start of the genome.
   *
   * @param g genome from organism being evaluated.
   * @return phenotype vector that is calculated from 'g'.
   */
  phenotype_t MultiPathExploration(const genome_t& g);

  /**
   * Multi-valley Crossing function:
   *
   * Solutions are pressured to cross valleys of different widths at each gene.
   * We use a peaks vector that supplied to calculate the penalty at each valley.
   * We use the difference between the gene value and penalty value to determine
   * the reduction when calculating the trait value.
   *
   * phenotype[i] = peaks[ floor(genome[i]) ] -  (genome[i] - peaks[ floor(genome[i]) ])
   *
   * @param g genome from organism being evaluated.
   * @param peaks vector of peaks for each floored gene value
   * @param dips_start gene value where peaks begin
   * @param dips_end gene value where dips end
   *
   * @return phenotype vector that is calculated from 'g'.
   */
  phenotype_t MultiValleyCrossing(
    const genome_t& g,
    const phenotype_t& peaks,
    double dips_start,
    double dips_end
  );

  ///< Functions that deal with interpretation of phenotype vectors

  /**
   * Optimized Vector:
   *
   * Will return a boolean vector that lists what genes are optimized relative to the target vector.
   * A gene is optimized if it meets the target value accuracy % requirement.
   *
   * @param g genome from organism being evaluated.
   * @param acc This value is the accuracy % needed to be considered optimized
   *
   * @return boolean vector that displays optimized traits.
  */
  emp::vector<bool> OptimizedVector(const genome_t& g, const double acc);

};

///< diagnostic problem implementations

Diagnostic::phenotype_t Diagnostic::ExploitationRate(const genome_t& g) {
  emp_assert(g.size() > 0);
  // Translate genotype directly into phenotype
  return phenotype_t(g);
}

Diagnostic::phenotype_t Diagnostic::OrderedExploitation(const genome_t & g) {
  emp_assert(g.size() > 0);
  emp_assert(cred_set);

  // find where descending order breaks
  const auto it = std::is_sorted_until(
    g.begin(),
    g.end(),
    std::greater<>()
  );

  // if sorted, return same vector; otherwise, fill out phenotype appropriately
  if(it == g.end()) {
    return phenotype_t(g);
  } else {
    // initialize phenotype to the same size
    phenotype_t phenotype(g.size());

    // calculate cutoff point where descending order is broken
    const size_t cutoff = std::distance(g.begin(), it);

    // everything up to unsorted
    for(size_t i = 0; i < cutoff; ++i) {
      phenotype[i] = g[i];
    }

    // everything after unsorted
    for(size_t i = cutoff; i < phenotype.size(); ++i) {
      phenotype[i] = max_cred;
    }

    return phenotype;
  }
}

Diagnostic::phenotype_t Diagnostic::ContradictoryObjectives(const genome_t& g) {
  emp_assert(g.size() > 0);
  emp_assert(cred_set);

  // intialize phenotype vector to max_cred in all positions
  phenotype_t phenotype(g.size(), max_cred);

  // find max value position
  const size_t max_gene = std::distance(
    g.begin(),
    std::max_element(g.begin(), g.end())
  );
  // set max position to gene value
  phenotype[max_gene] = g[max_gene];

  return phenotype;
}

Diagnostic::phenotype_t Diagnostic::MultiPathExploration(const genome_t& g) {
  emp_assert(0 < g.size());
  emp_assert(cred_set);

  // Intialize vector with size g, set each position to max_cred
  phenotype_t phenotype(g.size(), max_cred);

  // find max value position
  const auto opti_it = std::max_element(g.begin(), g.end());
  const size_t opti = std::distance(g.begin(), opti_it);

  // find where order breaks
  const size_t sort = std::distance(
    g.begin(),
    std::is_sorted_until(opti_it, g.end(), std::greater<>())
  );

  // middle of optimal value till order broken
  for(size_t i = opti; i < sort; ++i) {
    phenotype[i] = g[i];
  }

  return phenotype;
}

Diagnostic::phenotype_t Diagnostic::MultiValleyCrossing(
  const phenotype_t& p,
  const phenotype_t& peaks,
  double dips_start,
  double dips_end
) {
  // quick checks
  emp_assert(p.size() > 0);
  emp_assert(cred_set);
  emp_assert(peaks.size() > 0);

  // intialize vector with size p
  phenotype_t phenotype(p.size());

  for(size_t i = 0; i < p.size(); ++i) {
    phenotype[i] = (p[i] <= dips_start || p[i] >= dips_end) ? p[i] : 2.0 * peaks[static_cast<size_t>(p[i])] - p[i];

    // if (p[i] <= dips_start || p[i] >= dips_end) {
    //   phenotype[i] = p[i];
    // } else {
    //   phenotype[i] = 2.0 * peaks[static_cast<size_t>(p[i])] - p[i];
    // }
  }

  return phenotype;
}

///< phenotype vector interpretation implementations

emp::vector<bool> Diagnostic::OptimizedVector(
  const genome_t& g,
  const double acc
) {
  // quick checks
  emp_assert(g.size() > 0);
  emp_assert(target > 0);
  emp_assert(acc > 0);
  emp_assert(acc <= 1.0);

  // initialize optimized vector with all false
  emp::vector<bool> optimized(g.size(), false);
  const double threshold = acc * target;

  // iterate through genome and check optimality
  for(size_t i = 0; i < g.size(); ++i)
  {
    if (g[i] >= threshold) {
      optimized[i] = true;
    }
  }

  return optimized;
}

}

#endif