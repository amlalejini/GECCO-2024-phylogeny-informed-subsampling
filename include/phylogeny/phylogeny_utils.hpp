#pragma once

#include <map>
#include <set>
#include <unordered_set>
#include <optional>
#include <functional>
#include <utility>
#include <deque>

#include "emp/base/Ptr.hpp"
#include "emp/Evolve/Systematics.hpp"
#include "emp/base/vector.hpp"

// TODO - Data structure for tracking phenotypes
// TODO - Subclass of systematics that makes it easy to search for relative
// TODO - make search more efficient!

namespace phylo {

// NOTE - might need to make this more generic!
struct phenotype_info {

  using phen_t = emp::vector<double>;
  using has_phen_t = std::true_type;
  using has_mutations_t = std::false_type;
  using has_fitness_t = std::true_type;

  // Phenotype is set of test case scores
  // Need to set an eval vector

  emp::DataNode<
    double,
    emp::data::Current,
    emp::data::Range
  > fitness;            ///< Tracks fitness range for this taxon.

  phen_t phenotype;     ///< Tracks phenotype associated with this taxon.
  emp::vector<bool> traits_evaluated; ///< Tracks whether a particular trait has been evaluated

  // TODO - add struct that tracks each trait?
  // emp::vector<bool> traits_estimated;
  // emp::vector<double> trait_estimates;
  // emp::vector<size_t> estimate_dists;

  const phen_t& GetPhenotype() const {
    return phenotype;
  }

  double GetFitness() const {
    return fitness.GetMean();
  }

  const emp::vector<bool>& GetTraitsEvaluated() const {
    return traits_evaluated;
  }

  bool TraitEvaluated(size_t i) const {
    emp_assert(i < traits_evaluated.size());
    return traits_evaluated[i];
  }

  void RecordFitness(double fit) {
    fitness.Add(fit);
  }

  void RecordPhenotype(
    const phen_t& phen,
    const emp::vector<bool>& eval
  ) {
    phenotype = phen;
    traits_evaluated = eval;
  }

};

/// Return the nearest ancestor with the given trait evaluated.
template<typename TAXON>
std::optional<emp::Ptr<TAXON>> NearestAncestorWithTraitEval(
  emp::Ptr<TAXON> tax,
  size_t trait_id,
  size_t max_dist = (size_t)-1
) {
  emp_assert(tax != nullptr);
  size_t dist = 0;
  emp::Ptr<TAXON> cur_tax = tax;
  bool found_ancestor = false;

  // Keep searching unless we hit the root of the tree
  while (cur_tax != nullptr && dist <= max_dist) {
    // Localize relevant taxon information
    const auto& taxon_info = cur_tax->GetData();
    const auto& traits_evaluated = taxon_info.GetTraitsEvaluated();
    emp_assert(trait_id < traits_evaluated.size());
    // Has the focal trait been evaluted in this taxon?
    found_ancestor = traits_evaluated[trait_id];
    // If we found an ancestor with the requested trait evaluated OR we've
    // gone too far in the tree, break out of our search.
    if (found_ancestor) {
      break;
    }
    // Move up the tree, increase current search distance.
    cur_tax = cur_tax->GetParent();
    ++dist;
  }
  // If we found an ancestor, return it; otherwise, return empty optional.
  return (found_ancestor) ? std::optional<emp::Ptr<TAXON>>{cur_tax} : std::nullopt;
}

/// Return nearest relative with the given trait evaluated.
/// Breadth first traversal.
template<typename TAXON>
std::optional<emp::Ptr<TAXON>> NearestRelativeWithTraitEval(
  emp::Ptr<TAXON> tax,
  size_t trait_id,
  size_t max_dist = (size_t)-1
) {
  emp_assert(tax != nullptr);
  // Track discovered taxa. Start with current taxon.
  std::unordered_set<size_t> discovered_taxa({tax->GetID()});
  std::deque<
    std::pair<emp::Ptr<TAXON>, size_t>
  > search_queue;
  search_queue.emplace_back(tax, 0);

  while (!search_queue.empty()) {
    // Grab current taxon and its distance from search source.
    emp::Ptr<TAXON> cur_tax = search_queue.front().first;
    size_t cur_dist = search_queue.front().second;
    // Pop current taxon from front of search queue
    search_queue.pop_front();
    // Localize relevant taxon info
    const auto& traits_evaluated = cur_tax->GetData().GetTraitsEvaluated();
    emp_assert(trait_id < traits_evaluated.size(), trait_id, traits_evaluated.size());

    // Check if evaluated. If so, return current taxon.
    if (traits_evaluated[trait_id]) {
      return std::optional<emp::Ptr<TAXON>>{cur_tax};
    }

    // Add next set of relatives if not too far away.
    // Stop searching further if cur_dist >= max_dist
    if (cur_dist + 1 > max_dist) {
      continue;
    }

    // Add ancestors not already searched and not too far away
    emp::Ptr<TAXON> ancestor_tax = cur_tax->GetParent();
    if (ancestor_tax != nullptr && !emp::Has(discovered_taxa, ancestor_tax->GetID())) {
      discovered_taxa.emplace(ancestor_tax->GetID());
      search_queue.emplace_back(std::make_pair(ancestor_tax, cur_dist + 1));
    }

    // Add any descendants not already searched and not too far away
    std::set<emp::Ptr<TAXON>> descendants(cur_tax->GetOffspring());
    for (emp::Ptr<TAXON> descendant_tax : descendants) {
      emp_assert(descendant_tax != nullptr);
      if (!emp::Has(discovered_taxa, descendant_tax->GetID())) {
        discovered_taxa.emplace(descendant_tax->GetID());
        search_queue.emplace_back(std::make_pair(descendant_tax, cur_dist + 1));
      }
    }

  }
  return std::nullopt;
}

} // end namespace phylo