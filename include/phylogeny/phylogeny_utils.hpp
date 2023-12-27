#pragma once

#include <map>
#include <set>
#include <unordered_set>
#include <optional>
#include <functional>
#include <utility>
#include <deque>
#include <algorithm>

#include "emp/base/Ptr.hpp"
#include "emp/Evolve/Systematics.hpp"
#include "emp/base/vector.hpp"

// TODO - Data structure for tracking phenotypes
// TODO - Subclass of systematics that makes it easy to search for relative
// TODO - make search more efficient!

// * what population members are having the phylogenetic approximation applied,
// *what taxon is being used for the approximation, and
// *how far they are apart phylogeneticallly.

// TODO - after a few feature additions, these search functions need a redesign!

namespace phylo {

struct TraitEstInfo {
  bool estimated = false;         ///< Whether or not this trait has been searched from this taxon
  bool estimate_success = false;  ///< Whether or not the trait was successfully estimated
  double estimated_score = 0.0;   ///< Most recently estimated score
  size_t estimation_dist = 0;     ///< Distance of estimate
  size_t source_taxon_id = 0;     ///< Source taxon ID for the estimate

  TraitEstInfo(
    bool _estimated=false,
    bool _estimate_success=false,
    double _estimated_score=0.0,
    size_t _estimation_dist=0,
    size_t _source_taxon_id=0
  ) :
    estimated(_estimated),
    estimate_success(_estimate_success),
    estimated_score(_estimated_score),
    estimation_dist(_estimation_dist),
    source_taxon_id(_source_taxon_id)
  { ; }

  TraitEstInfo(const TraitEstInfo&) = default;
  TraitEstInfo(TraitEstInfo&&) = default;

  TraitEstInfo& operator=(const TraitEstInfo&) = default;
  TraitEstInfo& operator=(TraitEstInfo&&) = default;

  void SetEst(bool success, double score, size_t dist, size_t source_id) {
    estimated = true;
    estimate_success = success;
    estimated_score = score;
    estimation_dist = dist;
    source_taxon_id = source_id;
  }
};

// NOTE - might need to make this more generic!
struct taxon_info {

  using phen_t = emp::vector<double>;
  using has_phen_t = std::true_type;
  using has_mutations_t = std::true_type;
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
  bool phen_recorded=false;

  // TODO - implement rough mutation tracking!
  std::unordered_map<std::string, int> mut_counts;

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
    // If phenotype not recorded on this taxon, update phenotype.
    // Otherwise, merge incoming phenotype with previous phenotype.
    if (!phen_recorded) {
      phenotype = phen;
      traits_evaluated = eval;
      trait_est_info.clear();
      trait_est_info.resize(phenotype.size(), {false, false, 0, 0, 0});
      phen_recorded = true;
    } else {
      emp_assert(phen.size() == phenotype.size());
      emp_assert(eval.size() == traits_evaluated.size());
      for (size_t trait_id = 0; trait_id < phen.size(); ++trait_id) {
        // If this trait was already evaluated, don't update.
        if (traits_evaluated[trait_id]) continue;
        phenotype[trait_id] = phen[trait_id];
        traits_evaluated[trait_id] = eval[trait_id];
      }
      // Reset the trait estimation info
      emp_assert(trait_est_info.size() == phen.size());
      std::fill(
        trait_est_info.begin(),
        trait_est_info.end(),
        TraitEstInfo(false, false, 0, 0, 0)
      );
    }
  }

};

namespace impl {

void FilterEvaluated(size_t min_pool_size, emp::vector<size_t>& pool, const emp::vector<bool>& evaluated) {
  // Loop over pool (backwards), removing items in pool that have been evaluated
  for (int i = pool.size() - 1; i <= 0 && pool.size() > min_pool_size; --i) {
    emp_assert(i < (int)pool.size(), i, pool.size());
    const size_t trait_id = pool[(size_t)i];
    emp_assert(trait_id < evaluated.size());
    if (evaluated[trait_id]) {
      // If trait has been evaluated on this taxon:
      // - remove it from available
      // - add it to excluded
      pool.pop_back();
    }
  }
  emp_assert(pool.size() >= min_pool_size);
}

template<typename TAXON>
emp::vector<size_t> PhyloInformedSample_Ancestors(
  emp::Random& rnd,
  size_t sample_size,
  const emp::vector<size_t>& trait_ids,
  emp::Ptr<TAXON> tax,
  size_t max_dist = (size_t)-1
) {
  emp_assert(trait_ids.size() > 0);
  emp_assert(sample_size <= trait_ids.size(), "Number of choices must be at least as large as requested sample size");
  emp_assert(tax != nullptr);

  // Construct set of available choices
  emp::vector<size_t> available(trait_ids);
  emp::Shuffle(rnd, available);
  // Maintain vector of excluded traits (in order of their exclusion)
  emp::vector<size_t> excluded;

  // Keep winnowing down available choices until left with sample size
  size_t dist = 0;
  emp::Ptr<TAXON> cur_tax = tax;
  while (cur_tax != nullptr && dist <= max_dist) {
    // Localize relevant taxon information
    const auto& taxon_info = cur_tax->GetData();
    const auto& traits_evaluated = taxon_info.GetTraitsEvaluated();
    // If this taxon hasn't been evaluated yet, skip over
    if (traits_evaluated.size() == 0) {
      cur_tax = cur_tax->GetParent();
      ++dist;
      continue;
    }

    FilterEvaluated(sample_size, available, traits_evaluated);

    emp_assert(available.size() >= sample_size);
    if (available.size() == sample_size) {
      break;
    }

    // Move up to parent
    cur_tax = cur_tax->GetParent();
    ++dist;
  }
  // What happens if we hit max depth or we hit root & available still has too many things?
  // That means that *everything* in available hasn't been evaluated in an ancestor
  // => Resize down to sample size (random selection because we shuffled earlier)
  if (available.size() > sample_size) {
    available.resize(sample_size);
  }

  emp_assert(available.size() == sample_size);
  return available;
}

template<typename TAXON>
emp::vector<size_t> PhyloInformedSample_Relatives(
  emp::Random& rnd,
  size_t sample_size,
  const emp::vector<size_t>& trait_ids,
  emp::Ptr<TAXON> tax,
  size_t max_dist = (size_t)-1
) {
  emp_assert(trait_ids.size() > 0);
  emp_assert(sample_size <= trait_ids.size(), "Number of choices must be at least as large as requested sample size");
  emp_assert(tax != nullptr);

  // Construct set of available choices
  emp::vector<size_t> available(trait_ids);
  // Shuffle to eliminate bias
  emp::Shuffle(rnd, available);

  // Perform a breadth first search in phylogeny, winnowing down the
  // available choices until left with appropriate sample size

  // Track discovered taxa. Start with current taxon.
  std::unordered_set<size_t> discovered_taxa({tax->GetID()});
  std::deque<
    std::pair<emp::Ptr<TAXON>, size_t>
  > search_queue;
  search_queue.emplace_back(tax, 0);
  while (!search_queue.empty()) {
    // Grab current taxon and its distance from search source
    emp::Ptr<TAXON> cur_tax = search_queue.front().first;
    const size_t cur_dist = search_queue.front().second;
    // Pop current taxon from front of search queuej
    search_queue.pop_front();
    // Localize relevant taxon info
    const auto& taxon_info = cur_tax->GetData();
    const auto& traits_evaluated = taxon_info.GetTraitsEvaluated();

    // If this taxon has been evaluated, filter available choices
    if (traits_evaluated.size() > 0) {
      FilterEvaluated(sample_size, available, traits_evaluated);
      emp_assert(available.size() >= sample_size);
      if (available.size() == sample_size) {
        break;
      }
    }

    // Add next set of relatives to search queue if they aren't too far away
    if (cur_dist + 1 > max_dist) {
      continue;
    }

    // Add ancestors not already search and not too far away
    emp::Ptr<TAXON> ancestor_tax = cur_tax->GetParent();
    if (ancestor_tax != nullptr && !emp::Has(discovered_taxa, ancestor_tax->GetID())) {
      discovered_taxa.emplace(ancestor_tax->GetID());
      search_queue.emplace_back(
        std::make_pair(ancestor_tax, cur_dist + 1)
      );
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

  // What happens if we hit max depth or we hit root & available still has too many things?
  // That means that *everything* in available hasn't been evaluated in a relative
  // => Resize down to sample size (random selection because we shuffled earlier)
  if (available.size() > sample_size) {
    available.resize(sample_size);
  }
  emp_assert(available.size() == sample_size);
  return available;
}

}

// NOTE: The phylogeny annotations could be tweaked to further optimize the sampling
//       process. The current implementation is layered on top of pre-existing annotations
//       to minimize required changes to the tracking system.
//       If we wanted to refocus experiments on *just* phylogeny-based sampling,
//       it would be worthwhile to re-implement/re-factor with that in mind.
template<typename TAXON>
emp::vector<size_t> PhyloInformedSample(
  emp::Random& rnd,
  size_t sample_size,
  const emp::vector<size_t>& trait_ids,
  emp::Ptr<TAXON> tax,
  bool ancestors_only = true,
  size_t max_dist = (size_t)-1
) {
  return (ancestors_only) ?
    impl::PhyloInformedSample_Ancestors(rnd, sample_size, trait_ids, tax, max_dist) :
    impl::PhyloInformedSample_Relatives(rnd, sample_size, trait_ids, tax, max_dist);
}

/// Return the nearest ancestor with the given trait evaluated.
template<typename TAXON>
TraitEstInfo& NearestAncestorWithTraitEval(
  emp::Ptr<TAXON> tax,
  size_t trait_id,
  size_t max_dist = (size_t)-1,
  double est_fail_trait_val = 0
) {
  emp_assert(tax != nullptr);
  size_t dist = 0;
  emp::Ptr<TAXON> cur_tax = tax;
  TraitEstInfo& focal_est_info = tax->GetData().GetTraitEstimationInfo(trait_id);
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

  // If we found an ancestor, mark search as successful.
  if (found_ancestor) {
    focal_est_info.SetEst(
      true,
      cur_tax->GetData().GetTraitScore(trait_id),
      dist,
      cur_tax->GetID()
    );
  } else {
    focal_est_info.SetEst(
      false,
      est_fail_trait_val,
      0,
      tax->GetID()
    );
  }

  return focal_est_info;
}

/// Return the nearest ancestor with the given trait evaluated.
template<typename TAXON>
TraitEstInfo& NearestAncestorWithTraitEvalOpt(
  emp::Ptr<TAXON> tax,
  size_t trait_id,
  size_t max_dist = (size_t)-1,
  double est_fail_trait_val = 0.0
) {
  emp_assert(tax != nullptr);
  size_t dist = 0;
  emp::Ptr<TAXON> cur_tax = tax;
  bool found_ancestor = false;
  double est_score = 0.0;
  size_t est_id = 0;

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
      est_id = cur_tax->GetID();
      break;
    }

    // Has the focal trait been estimated in this taxon before?
    const auto& trait_est_info = taxon_info.GetTraitEstimationInfo(trait_id);
    // if (trait_est_info.estimated) {
    if (trait_est_info.estimate_success) {
      // Only mark ancestor as found if estimation falls within distance constraint
      found_ancestor = ((trait_est_info.estimation_dist + dist) <= max_dist);
      est_score = trait_est_info.estimated_score;
      dist = trait_est_info.estimation_dist + dist;
      est_id = trait_est_info.source_taxon_id;
      break;
    }

    // Move up the tree, increase current search distance.
    cur_tax = cur_tax->GetParent();
    ++dist;
  }


  // If we found an ancestor with the trait, note it in the given taxon.
  auto& focal_taxon_info = tax->GetData();
  auto& focal_est_info = focal_taxon_info.GetTraitEstimationInfo(trait_id);
  if (found_ancestor) {
    focal_est_info.SetEst(
      true,
      est_score,
      dist,
      est_id
    );
  } else {
    focal_est_info.SetEst(
      false,
      est_fail_trait_val,
      0,
      tax->GetID()
    );
  }

  // If we found an ancestor, return it; otherwise, return empty optional.
  return focal_est_info;
}

/// Return nearest relative with the given trait evaluated.
/// Breadth first traversal.
template<typename TAXON>
TraitEstInfo& NearestRelativeWithTraitEval(
  emp::Ptr<TAXON> tax,
  size_t trait_id,
  size_t max_dist = (size_t)-1,
  double est_fail_trait_val = 0.0
) {
  emp_assert(tax != nullptr);
  // Track discovered taxa. Start with current taxon.
  std::unordered_set<size_t> discovered_taxa({tax->GetID()});
  std::deque<
    std::pair<emp::Ptr<TAXON>, size_t>
  > search_queue;
  search_queue.emplace_back(tax, 0);

  TraitEstInfo& focal_est_info = tax->GetData().GetTraitEstimationInfo(trait_id);

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
      focal_est_info.SetEst(
        true,
        cur_tax->GetData().GetTraitScore(trait_id),
        cur_dist,
        cur_tax->GetID()
      );
      return focal_est_info;
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

  focal_est_info.SetEst(
    false,
    est_fail_trait_val,
    0,
    tax->GetID()
  );
  return focal_est_info;
}


/// Return nearest relative with the given trait evaluated.
/// Breadth first traversal.
template<typename TAXON>
TraitEstInfo& NearestRelativeWithTraitEvalOpt(
  emp::Ptr<TAXON> tax,
  size_t trait_id,
  size_t max_dist = (size_t)-1,
  double est_fail_trait_val = 0.0
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
  size_t est_id = 0;

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
      est_id = cur_tax->GetID();
      break;
    }

    // If relative has been estimated & estimation is too far away, skip.
    auto& cur_tax_info = cur_tax->GetData();
    auto& cur_tax_est = cur_tax_info.GetTraitEstimationInfo(trait_id);

    // if (cur_tax_est.estimated && cur_tax_est.estimation_dist >= max_dist) {
    if (cur_tax_est.estimate_success && cur_tax_est.estimation_dist >= max_dist) {
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

  if (found_relative) {
    focal_est_info.SetEst(
      true,
      est_score,
      est_dist,
      est_id
    );
  } else {
    focal_est_info.SetEst(
      false,
      est_fail_trait_val,
      0,
      tax->GetID()
    );
  }

  return focal_est_info;
}

} // end namespace phylo
