/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/

#ifndef DIAGNOSTICS_ORG_HPP
#define DIAGNOSTICS_ORG_HPP

#include <algorithm>
#include <numeric>
#include "emp/base/vector.hpp"

namespace diag {

// TODO - move functions that should be protected into protected scope
// TODO - clean up function interface
class Org {
public:
  using genome_t = emp::vector<double>;     // genome vector type
  using phenotype_t = emp::vector<double>;  // score vector type

  static constexpr double START_DB = 0.0;    // Initial gene values
  static constexpr size_t START_ST = 0;

protected:
  genome_t genome; // organism genome vector
  phenotype_t phenotype;        // organism score vector
  bool evaluated = false; // has this organism been evaluated yet?


  bool scored = false;  // score vector set?

  emp::vector<bool> optimal_genes; // organims gene optimal vector
  bool optimality_calculated = false; // gene optimal vector calculated?

  size_t optimal_gene_count = 0;     // optimal gene count
  bool optimal_genes_counted = false; // gene optimal vector counted?

  double aggregate_score = 0.0; // aggregate score
  bool aggregated = false; // aggregate calculate?

  size_t streak = 0;     // streak count
  bool streaked = false; // streak calculated?

  size_t M = 0;         // Number of genes in genome

  size_t start_pos;   // starting position
  bool start = false; // starting position located?

  bool clone = false; // Are we a clone?

  // set the aggregated phenotype (called from world.h or inherited from parent)
  void SetAggregate(double agg) {
    emp_assert(!aggregated);
    emp_assert(M > 0);
    aggregated = true;
    aggregate_score = agg;
  }

  // set the optimal gene count (called from world.h or inherited from parent)
  void SetOptimalGeneCount(size_t count) {
    emp_assert(optimality_calculated);
    emp_assert(!optimal_genes_counted);
    emp_assert(M > 0);
    optimal_genes_counted = true;
    optimal_gene_count = count;
  }

  // set the starting position
  void SetStart(size_t pos) {
    emp_assert(!start);
    emp_assert(M > 0);
    start = true;
    start_pos = pos;
  }

  // set the starting position
  void SetStreak(size_t pos) {
    emp_assert(!streaked);
    emp_assert(0 < M);
    streaked = true;
    streak = pos;
  }

public:

  Org(size_t num_genes) {
    emp_assert(genome.size() == 0);
    emp_assert(M == 0);
    M = num_genes;
    start_pos = M;
    streak = M;
    genome.resize(M, Org::START_DB);
    optimal_genes.resize(M, false);
  }

  // every org after starting generation
  Org(const genome_t & _g) {
    emp_assert(genome.size() == 0);
    emp_assert(M == 0);
    M = _g.size();
    start_pos = _g.size();
    streak = _g.size();
    genome.resize(M, 0.0);
    std::copy(_g.begin(), _g.end(), genome.begin());
    optimal_genes.resize(M, false);
  }

  Org(const Org&) = default;
  Org(Org&&) = default;
  ~Org() { ; }
  Org& operator=(const Org&) = default;
  Org& operator=(Org&&) = default;

  // --- accessors ---
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
  phenotype_t & GetPhenotype() {
    emp_assert(evaluated);
    return phenotype;
  }

  const emp::vector<bool>& GetOptimalGenes() const {
    emp_assert(optimality_calculated);
    return optimal_genes;
  }
  emp::vector<bool>& GetOptimalGenes() {
    emp_assert(optimality_calculated);
    return optimal_genes;
  }

  double GetAggregateScore() const {
    emp_assert(aggregated);
    return aggregate_score;
  }

  bool IsClone() const {
    return clone;
  }

  size_t GetOptimalGeneCount() const {
    emp_assert(optimal_genes_counted);
    return optimal_gene_count;
  }

  size_t GetGeneCount() const {
    return M;
  }

  // get start position
  size_t GetStart() const {
    emp_assert(start_pos != M);
    return start_pos;
  }

  // get streak count
  size_t GetStreak() const {
    emp_assert(streak != M);
    return streak;
  }

  // Are we optimized at this objective?
  bool OptimizedAt(size_t obj) const;
  // get evaluated bool
  bool GetEvaluated() const {
    return evaluated;
  }
  // get optimal bool
  bool IsOptimalCalculated() const {
    return optimality_calculated;
  }
  // get aggregated bool
  bool GetAggregated() const {
    return aggregated;
  }
  // get counted bool
  bool IsOptimalGenesCounted() const {
    return optimal_genes_counted;
  }
  // get streak bool
  bool IsStreaked() const {
    return streaked;
  }
  // get max trait
  double GetMaxTrait() const {
    emp_assert(start);
    emp_assert(M > 0);
    emp_assert(phenotype.size() == M);
    return phenotype[start_pos];
  }

  // set phenotype vector (recieved from problem.h in world.h or inherited from parent)
  void SetPhenotype(const phenotype_t& phen) {
    // make sure that phenotype vector hasn't been set before.
    emp_assert(!evaluated);
    emp_assert(phen.size() == M);
    emp_assert(phenotype.size() == 0);
    emp_assert(M > 0);
    evaluated = true;
    phenotype = phen;
    AggregateScore();
    StartPosition();
  }

  // Set the optimal gene vector (recieved from problem.h in world.h or inherited from parent)
  void SetOptimal(const emp::vector<bool>& opt) {
    // make sure that optimal gene vector hasn't been set before.
    emp_assert(!optimality_calculated);
    emp_assert(opt.size() == M);
    emp_assert(optimal_genes.size() == 0);
    emp_assert(M > 0);
    optimality_calculated = true;
    optimal_genes = opt;
  }

  /**
   * Aggregate Score function:
   *
   * First we check to see if the phenotype has been calculated.
   * If yes, aggregated=True and throw expection.
   * Else, aggregated=False and we calculate it and return aggregate.
   *
   * @return aggregate
  */
  double AggregateScore();

  /**
   * Count Optimized function:
   *
   * First we check to see if the count has been calculated.
   * If yes, counted=True and throw expection.
   * Else, countedd=False and we calculate it and return count.
   *
   * @return count
  */
  size_t CountOptimizedGenes();

  /**
   * Find Starting Position
   *
   * Find the staring position in the phenotype vector.
   */
  size_t StartPosition();


  /**
   * Find Maximum Streak
   *
   * Find the biggest streak.
   */
  size_t CalcStreak();

  /**
   * Reset function:
   *
   * Will reset all variables in an organims when birth occurs.
   * Function executes if offspring is not a clone
  */
  void Reset();

  /**
   * Clone function:
   *
   * Will pass all info from parent to offspring solution.
   * Function executes if offpsring is an clone
   *
   * @param s phenotype vector recived
   * @param o optimal gene vector recieved
   * @param c optimal gene count recieved
   * @param a aggregate phenotype recieved
   *
  */
  void Clone(
    const phenotype_t& s,
    const emp::vector<bool>& o,
    size_t c,
    double a,
    size_t st,
    size_t sr
  );

};

bool Org::OptimizedAt(size_t obj) const {
  emp_assert(obj < M);
  emp_assert(optimal_genes.size() > 0);
  emp_assert(M == optimal_genes.size());

  return optimal_genes[obj];
}

double Org::AggregateScore() {
  emp_assert(!aggregated);
  emp_assert(M > 0);
  emp_assert(phenotype.size() == M, phenotype.size());

  // calculate the aggregate phenotype and set it
  const double agg = std::accumulate(
    phenotype.begin(),
    phenotype.end(),
    Org::START_DB
  );
  // Update internal aggregate score value
  SetAggregate(agg);
  return aggregate_score;
}

size_t Org::CountOptimizedGenes() {
  emp_assert(!optimal_genes_counted);
  emp_assert(M > 0);
  emp_assert(optimality_calculated);
  emp_assert(optimal_genes.size() == M, optimal_genes.size());

  // calculate total optimal genes and set it
  SetOptimalGeneCount(
    std::accumulate(
      optimal_genes.begin(),
      optimal_genes.end(),
      Org::START_ST
    )
  );

  return optimal_gene_count;
}

size_t Org::StartPosition() {
  emp_assert(!start);
  emp_assert(M > 0);
  emp_assert(phenotype.size() == M);

  // find max value position
  auto opti_it = std::max_element(phenotype.begin(), phenotype.end());
  SetStart(
    std::distance(phenotype.begin(), opti_it)
  );

  return start_pos;
}

size_t Org::CalcStreak() {
  emp_assert(!streaked);
  emp_assert(M > 0);
  emp_assert(phenotype.size() == M);

  // get longest streak
  size_t count = 0;
  size_t max_cnt = 0;
  for(auto& s : phenotype) {
    if(s > 0.0) {
      count++;
    } else {
      if(count > max_cnt) {
        max_cnt = count;
      }
      count = 0;
    }
  }

  // streak = max_cnt;
  SetStreak(max_cnt);
  return streak;
}

void Org::Reset() {
  // quick checks
  emp_assert(M > 0);
  emp_assert(genome.size() > 0);

  // reset phenotype vector
  phenotype.clear();
  evaluated = false;

  // reset optimal gene vector
  optimal_genes.clear();
  optimality_calculated = false;

  // reset optimal gene count stuff
  optimal_gene_count = 0;
  optimal_genes_counted = false;

  // reset aggregate phenotype stuff
  aggregate_score = 0.0;
  aggregated = false;

  // reset starting position info
  start_pos = genome.size();
  start = false;

  // reset clone var
  clone = false;
}

void Org::Clone(
  const phenotype_t& s,
  const emp::vector<bool>& o,
  size_t c,
  double a,
  size_t st,
  size_t sr
) {
  emp_assert(M > 0);
  emp_assert(genome.size() > 0);
  clone = true;
  // copy everything into offspring solution
  SetPhenotype(s);
  SetOptimal(o);
  SetOptimalGeneCount(c);
  SetAggregate(a);
  SetStart(st);
  SetStreak(sr);
}

} // namespace diag

#endif