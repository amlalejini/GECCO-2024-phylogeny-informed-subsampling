#pragma once

#include <algorithm>
#include <tuple>

#include "emp/base/vector.hpp"

#include "ProgSynthGenome.hpp"
#include "ProgSynthPhenotype.hpp"
#include "TestResult.hpp"

namespace psynth {

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