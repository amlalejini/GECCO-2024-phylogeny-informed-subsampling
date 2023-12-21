#pragma once

#include <unordered_map>
#include <string>

#include "../phylogeny/phylogeny_utils.hpp"

namespace psynth {

struct ProgSynthTaxonInfo : public phylo::phenotype_info {

  using phylo::phenotype_info::phen_t;
  using phylo::phenotype_info::has_phen_t;
  using phylo::phenotype_info::has_fitness_t;
  using has_mutations_t = std::true_type;

  std::unordered_map<std::string, int> mut_counts;

  // TODO

};


} // End psynth namespace