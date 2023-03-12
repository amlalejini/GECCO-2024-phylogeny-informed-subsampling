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

struct TraitEstInfo {
  size_t estimated = false;     ///< Whether or not this trait has an estimated score
  double estimated_score = 0.0; ///< Most recently estimated score
  size_t estimation_dist = 0;   ///< Distance of estimate

  TraitEstInfo(
    bool _estimated=false,
    double _estimated_score=0.0,
    size_t _estimation_dist=0
  ) :
    estimated(_estimated),
    estimated_score(_estimated_score),
    estimation_dist(_estimation_dist)
  { ; }

  TraitEstInfo(const TraitEstInfo&) = default;
  TraitEstInfo(TraitEstInfo&&) = default;

  TraitEstInfo& operator=(const TraitEstInfo&) = default;
  TraitEstInfo& operator=(TraitEstInfo&&) = default;
};

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

  phen_t phenotype;                   ///< Tracks phenotype associated with this taxon.
  emp::vector<bool> traits_evaluated; ///< Tracks whether a particular trait has been evaluated
  emp::vector<TraitEstInfo> trait_est_info;  ///< Tracks information about trait estimation.

  const phen_t& GetPhenotype() const {
    return phenotype;
  }

  double GetTraitScore(size_t i) const {
    return phenotype[i];
  }

  double GetFitness() const {
    return fitness.GetMean();
  }

  const emp::vector<bool>& GetTraitsEvaluated() const {
    return traits_evaluated;
  }

  const emp::vector<TraitEstInfo>& GetTraitsEstimationInfo() const {
    return trait_est_info;
  }

  TraitEstInfo& GetTraitEstimationInfo(size_t i) {
    emp_assert(i < trait_est_info.size());
    return trait_est_info[i];
  }

  const TraitEstInfo& GetTraitEstimationInfo(size_t i) const {
    emp_assert(i < trait_est_info.size());
    return trait_est_info[i];
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
    emp_assert(phen.size() == eval.size());
    phenotype = phen;
    traits_evaluated = eval;
    trait_est_info.clear();
    trait_est_info.resize(phenotype.size(), {false, 0, 0});
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

/// Return the nearest ancestor with the given trait evaluated.
template<typename TAXON>
std::optional<emp::Ptr<TAXON>> NearestAncestorWithTraitEvalOpt(
  emp::Ptr<TAXON> tax,
  size_t trait_id,
  size_t max_dist = (size_t)-1
) {
  emp_assert(tax != nullptr);
  size_t dist = 0;
  emp::Ptr<TAXON> cur_tax = tax;
  bool found_ancestor = false;
  double est_score = 0.0;

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
      est_score = taxon_info.GetTraitScore(trait_id);
      break;
    }

    // Has the focal trait been estimated in this taxon before?
    const auto& trait_est_info = taxon_info.GetTraitEstimationInfo(trait_id);
    if (trait_est_info.estimated) {
      // Only mark ancestor as found if estimation falls within distance constraint
      found_ancestor = ((trait_est_info.estimation_dist + dist) <= max_dist);
      est_score = trait_est_info.estimated_score;
      dist = trait_est_info.estimation_dist + dist;
      break;
    }

    // Move up the tree, increase current search distance.
    cur_tax = cur_tax->GetParent();
    ++dist;
  }

  // If we found an ancestor with the trait, note it in the given taxon.
  auto& focal_taxon_info = tax->GetData();
  auto& focal_est_info = focal_taxon_info.GetTraitEstimationInfo(trait_id);
  focal_est_info.estimated = true;
  focal_est_info.estimated_score = (dist <= max_dist) ? est_score : focal_taxon_info.GetTraitScore(trait_id);
  focal_est_info.estimation_dist = dist;

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


/// Return nearest relative with the given trait evaluated.
/// Breadth first traversal.
template<typename TAXON>
std::optional<emp::Ptr<TAXON>> NearestRelativeWithTraitEvalOpt(
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

  emp::Ptr<TAXON> cur_tax = nullptr;
  bool found_relative = false;
  double est_score = 0;
  size_t est_dist = 0;

  auto& focal_taxon_info = tax->GetData();
  auto& focal_est_info = focal_taxon_info.GetTraitEstimationInfo(trait_id);

  while (!search_queue.empty()) {
    // Grab current taxon and its distance from search source.
    cur_tax = search_queue.front().first;
    size_t cur_dist = search_queue.front().second;
    // Pop current taxon from front of search queue
    search_queue.pop_front();
    // Localize relevant taxon info
    const auto& traits_evaluated = cur_tax->GetData().GetTraitsEvaluated();
    emp_assert(trait_id < traits_evaluated.size(), trait_id, traits_evaluated.size());

    // Check if evaluated. If so, return current taxon.
    if (traits_evaluated[trait_id]) {
      found_relative = true;
      est_score = cur_tax->GetData().GetTraitScore(trait_id);
      est_dist = cur_dist;
      break;
    }

    // If relative has been estimated & estimation is too far away, skip.
    auto& cur_tax_info = cur_tax->GetData();
    auto& cur_tax_est = cur_tax_info.GetTraitEstimationInfo(trait_id);
    if (cur_tax_est.estimated && cur_tax_est.estimation_dist >= max_dist) {
      continue;
    }

    // Add next set of relatives if not too far away.
    // Stop searching further if cur_dist >= max_dist
    if (cur_dist + 1 > max_dist) {
      est_dist = cur_dist + 1;
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

  focal_est_info.estimated = true;
  focal_est_info.estimation_dist = est_dist;
  if (found_relative) {
    focal_est_info.estimated_score = est_score;
  } else {
    // If not relative found, just fill in with own score (probably 0).
    focal_est_info.estimated_score = focal_taxon_info.GetTraitScore(trait_id);
  }

  return (found_relative) ? std::optional<emp::Ptr<TAXON>>{cur_tax} : std::nullopt;
}

} // end namespace phylo