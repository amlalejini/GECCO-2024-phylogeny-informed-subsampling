#pragma once

#include <unordered_map>
#include <string>

#include "../phylogeny/phylogeny_utils.hpp"

namespace psynth {

struct ProgSynthTaxonInfo : public phylo::taxon_info {

  using phylo::taxon_info::phen_t;
  using phylo::taxon_info::has_phen_t;
  using phylo::taxon_info::has_fitness_t;
  using has_mutations_t = std::true_type;

  emp::vector<double> true_training_scores; ///< Used to store true training scores for a taxon when taking a phylogeny snapshot. Not used for search.
  double true_agg_score=0.0;
  bool true_training_scores_computed=false;

};


} // End psynth namespace