#pragma once

#include <functional>
#include <algorithm>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"

#include "BaseSelect.hpp"

namespace selection {

struct TournamentSelect : public BaseSelect {
public:
  using score_fun_t = std::function<double(void)>;
protected:
  emp::vector<score_fun_t>& score_funs;     ///< One function for each selection candidate (e.g., each member of the population)
  emp::Random& random;
  size_t tournament_size;

  // --- INTERNAL ---
  emp::vector<size_t> all_candidate_ids;

public:
  TournamentSelect(
    emp::vector<score_fun_t>& a_score_funs,
    emp::Random& a_random,
    size_t a_tournament_size=4
  ) :
    score_funs(a_score_funs),
    random(a_random),
    tournament_size(a_tournament_size)
  { ; }

  emp::vector<size_t>& operator()(size_t n) override;
  emp::vector<size_t>& operator()(
    size_t n,
    const emp::vector<size_t>& candidate_ids
  );

  size_t GetTournamentSize() const { return tournament_size; }
  void SetTournamentSize(size_t t) { tournament_size = t; }

};

emp::vector<size_t>& TournamentSelect::operator()(size_t n) {
  emp_assert(score_funs.size() > 0);
  const size_t num_candidates = score_funs.size();
  if (all_candidate_ids.size() != num_candidates) {
    all_candidate_ids.resize(num_candidates, 0);
    std::iota(
      all_candidate_ids.begin(),
      all_candidate_ids.end(),
      0
    );
  }
  return (*this)(n, all_candidate_ids);
}

emp::vector<size_t>& TournamentSelect::operator()(
  size_t n,
  const emp::vector<size_t>& candidate_ids
) {
  emp_assert(tournament_size > 0, "Tournament size must be greater than 0.", tournament_size);
  emp_assert(tournament_size <= score_funs.size());
  emp_assert(tournament_size <= candidate_ids.size());
  emp_assert(candidate_ids.size() <= score_funs.size());

  const size_t num_candidates = candidate_ids.size();
  selected.resize(n, 0);

  emp::vector<size_t> entries(tournament_size, 0);
  emp::vector<size_t> candidate_idxs(num_candidates, 0);
  std::iota(
    candidate_idxs.begin(),
    candidate_idxs.end(),
    0
  );

  for (size_t t = 0; t < n; ++t) {
    // Form a tournament
    emp::Shuffle(random, candidate_idxs);
    for (size_t i = 0; i < entries.size(); ++i) {
      entries[i] = candidate_ids[candidate_idxs[i]];
    }
    // pick a winner
    size_t winner_id = entries[0];
    double winner_fit = score_funs[entries[0]]();
    for (size_t i = 1; i < entries.size(); ++i) {
      const size_t entry_id = entries[i];
      const double entry_fit = score_funs[entry_id]();
      if (entry_fit > winner_fit) {
        winner_id = entry_id;
        winner_fit = entry_fit;
      }
    }
    // save the winner of tournament t
    selected[t] = winner_id;
  }
  return selected;
}

}