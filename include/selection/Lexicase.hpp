#pragma once

#include <functional>
#include <algorithm>

#include "emp/base/vector.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"

#include "BaseSelect.hpp"

namespace selection {

struct LexicaseSelect : public BaseSelect {
public:
  using score_fun_t = std::function<double(void)>;
protected:
  emp::vector< emp::vector<score_fun_t> >& score_fun_sets; ///< Per-candidate, per-function
  emp::Random& random;

  // --- INTERNAL ---
  emp::vector< emp::vector<double> > score_table;   ///< Used internally to track relevant scores for a given selection event

  emp::vector<size_t> score_ordering;               ///< Used internally to track function ordering. WARNING - Don't change values in this!
  emp::vector<size_t> candidate_idxs;               ///< Does not contain ids. Contains indices into score table.
  emp::vector<size_t> all_candidate_ids;
  emp::vector<size_t> all_fun_ids;
public:

  LexicaseSelect(
    emp::vector< emp::vector<score_fun_t> >& a_score_fun_sets,
    emp::Random& a_random
  ) :
    score_fun_sets(a_score_fun_sets),
    random(a_random)
  { ; }

  emp::vector<size_t>& operator()(size_t n) override;

  emp::vector<size_t>& operator()(
    size_t n,
    const emp::vector<size_t>& candidate_ids,
    const emp::vector<size_t>& fun_ids
  );

};

emp::vector<size_t>& LexicaseSelect::operator()(size_t n) {
  emp_assert(score_fun_sets.size() > 0);
  const size_t num_candidates = score_fun_sets.size();
  const size_t num_funs = score_fun_sets.back().size();
  // Use all candidates (fill out all candidate ids if doesn't match score fun set)
  if (all_candidate_ids.size() != num_candidates) {
    all_candidate_ids.resize(num_candidates, 0);
    std::iota(
      all_candidate_ids.begin(),
      all_candidate_ids.end(),
      0
    );
  }
  // Use all fitness functions (fill out function ids if doesn't match)
  if (all_fun_ids.size() != num_funs) {
    all_fun_ids.resize(num_funs, 0);
    std::iota(
      all_fun_ids.begin(),
      all_fun_ids.end(),
      0
    );
  }
  return (*this)(n, all_candidate_ids, all_fun_ids);
}

emp::vector<size_t>& LexicaseSelect::operator()(
  size_t n,
  const emp::vector<size_t>& candidate_ids,
  const emp::vector<size_t>& fun_ids
) {
  // TODO - test lexicase selection
  const size_t num_candidates = candidate_ids.size();
  const size_t fun_cnt = fun_ids.size();
  emp_assert(score_fun_sets.size() > 0);
  emp_assert(num_candidates > 0);
  emp_assert(num_candidates <= score_fun_sets.size());
  emp_assert(fun_cnt > 0);
  emp_assert(fun_cnt <= score_fun_sets.back().size());

  // Reset internal selected vector to size n
  selected.resize(n, 0);

  // Update the score table
  score_table.resize(fun_cnt);
  for (size_t fun_i = 0; fun_i < fun_cnt; ++fun_i) {
    score_table[fun_i].resize(num_candidates);
    const size_t fun_id = fun_ids[fun_i];

    for (size_t cand_i = 0; cand_i < num_candidates; ++cand_i) {
      emp_assert(fun_cnt <= score_fun_sets[cand_i].size());
      emp_assert(fun_id < score_fun_sets[cand_i].size());
      const size_t cand_id = candidate_ids[cand_i];
      score_table[fun_i][cand_i] = score_fun_sets[cand_id][fun_id]();
    }
  }

  // Update score ordering
  if (fun_cnt != score_ordering.size()) {
    score_ordering.resize(fun_cnt);
    std::iota(
      score_ordering.begin(),
      score_ordering.end(),
      0
    );
  }
  // Update valid candidate indices into score_table
  if (num_candidates != candidate_idxs.size()) {
    candidate_idxs.resize(num_candidates);
    std::iota(
      candidate_idxs.begin(),
      candidate_idxs.end(),
      0
    );
  }
  emp::vector<size_t> cur_pool, next_pool;
  for (size_t sel_i = 0; sel_i < n; ++sel_i) {
    // Randomize the score ordering
    emp::Shuffle(random, score_ordering);
    // Step through each score
    cur_pool = candidate_idxs;
    int depth = -1;
    // For each score, filter the population down to only the best performers.
    for (size_t score_id : score_ordering) {
      ++depth;
      double max_score = score_table[score_id][cur_pool[0]]; // Max score starts as first candidate's score on this function.
      next_pool.emplace_back(cur_pool[0]); // Seed the keeper pool with the first candidate.

      for (size_t i = 1; i < cur_pool.size(); ++i) {
        const size_t cand_idx = cur_pool[i];
        const double cur_score = score_table[score_id][cand_idx];
        if (cur_score > max_score) {
          max_score = cur_score;        // This is the new max score for this function
          next_pool.resize(1);          // Clear out candidates with former max score
          next_pool[0] = cand_idx;      // Add this candidate as only one with the new max
        } else if (cur_score == max_score) {
          next_pool.emplace_back(cand_idx);
        }
      }
      // Make next_pool into new cur_pool; make cur_pool allocated space for next_pool
      std::swap(cur_pool, next_pool);
      next_pool.resize(0);
      if (cur_pool.size() == 1) break; // Stop if we're down to just one candidate.
    }
    // Select a random survivor (all equal at this point)
    emp_assert(cur_pool.size() > 0);
    const size_t win_id = (cur_pool.size() == 1)
      ? cur_pool.back()
      : cur_pool[random.GetUInt(cur_pool.size())];
    emp_assert(win_id < candidate_ids.size());
    selected[sel_i] = candidate_ids[win_id]; // Transform win id (which is an index into score table) back into an actual candidate id
  }
  return selected;
}

} // End selection namespace