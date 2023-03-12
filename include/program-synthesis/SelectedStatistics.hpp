#pragma once

#include <utility>
#include <algorithm>

#include "emp/base/vector.hpp"
#include "emp/math/stats.hpp"

namespace psynth {

  struct SelectedStatistics {
    size_t num_unique_cand_selected;
    double entropy_cand_selected;
    size_t parents_num_tests_covered;
    emp::vector<bool> parent_test_coverage;

    void Reset() {
      parent_test_coverage.clear();
      num_unique_cand_selected = 0;
      entropy_cand_selected = 0;
      parents_num_tests_covered = 0;
    }

    template<typename WORLD_T>
    void Calculate(
      const emp::vector<size_t>& selected,
      WORLD_T& world
    ) {
      Reset();
      num_unique_cand_selected = emp::UniqueCount(selected);
      entropy_cand_selected = emp::ShannonEntropy(selected);
      size_t num_tests = world.GetNumTrainingCases();
      parent_test_coverage.resize(num_tests, false);
      for (size_t test_id = 0; test_id < num_tests; ++test_id) {
        for (size_t s_i = 0; (s_i < selected.size() && !parent_test_coverage[test_id]); ++s_i) {
          const size_t org_id = selected[s_i];
          auto& org = world.GetOrg(org_id);
          const bool pass_test = org.GetPhenotype().PassedTest(test_id);
          parent_test_coverage[test_id] = parent_test_coverage[test_id] || pass_test;
        }
      }
      parents_num_tests_covered = std::accumulate(
        parent_test_coverage.begin(),
        parent_test_coverage.end(),
        0
      );
    }
  };

}