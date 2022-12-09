/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/
#pragma once

#include <iostream>

#include "emp/Evolve/World.hpp"

#include "phylogeny/Phylogeny.hpp"
#include "selection/SelectionSchemes.hpp"

#include "DiagnosticsConfig.hpp"
#include "DiagnosticsOrg.hpp"
#include "problems/DiagnosticsProblems.hpp"

namespace diag {

class DiagnosticsWorld : public emp::World<DiagnosticsOrg> {
public:
  using base_t = emp::World<DiagnosticsOrg>;
  using org_t = DiagnosticsOrg;
  using genome_t = typename org_t::genome_t;
  using phenotype_t = typename org_t::phenotype_t;

  using config_t = DiagnosticsConfig;
  using selection_fun_t = std::function< emp::vector<size_t>&(size_t, const emp::vector<size_t>&, const emp::vector<size_t>&) >;

  struct Grouping {
    size_t group_id = 0;
    size_t group_size = 0;
    emp::vector<size_t> member_ids;

    void Resize(size_t size, size_t value=0) {
      group_size = size;
      member_ids.resize(group_size, value);
    }

    size_t GetSize() const { return member_ids.size(); }
  };

protected:
  const config_t& config;

  emp::Ptr<BaseDiagnostic> base_diagnostic=nullptr; ///< Base-layer diagnostic to use to translate genomes to phenotypes
  MultiValleyCrossingDiagnostic valley_diagnostic;  ///< If in use, it's layered on top of another diagnostic (base_diagnostic).

  // std::function<void(org_t&)> evaluate_org_fun;
  emp::Signal<void(size_t)> do_org_evaluation_sig;
  std::function<void(const genome_t&, phenotype_t&)> translate_genome_fun;

  // TODO - create a class/struct that manages all of this?
  size_t total_tests=0;
  emp::vector<double> org_aggregate_scores;
  // emp::vector< std::function<double(const org_t&)> > fit_fun_set; ///< One function for every possible test case.
  emp::vector< emp::vector< std::function<double(void)> > > fit_fun_set; ///< Per-organism, per-test
  emp::vector< std::function<double(void)> > agg_score_fun_set; ///< Per-organism, aggregate score

  emp::vector< emp::vector<double> > org_test_scores;   ///< Test scores for each organism
  emp::vector< emp::vector<bool> > org_test_evaluations; ///< Which test cases has each organism been evaluated on?

  emp::vector<size_t> possible_test_ids;
  emp::vector<size_t> possible_pop_ids;

  emp::vector<Grouping> test_groupings; ///< Groupings of tests (needs to be same size as org_groupings)
  emp::vector<Grouping> org_groupings;  ///< Groupings of organisms (needs to be same size as test_groupings)
  emp::vector<size_t> org_group_ids;    ///< Group ID for each organism

  std::function<void(void)> assign_test_groupings;
  std::function<void(void)> assign_org_groupings;

  std::function<double(size_t, size_t)> estimate_test_score;

  emp::Ptr<selection::BaseSelect> selector;
  emp::vector<size_t> selected_parent_ids;

  std::function<void(void)> run_selection_routine;
  selection_fun_t selection_fun;

  size_t total_test_evaluations = 0;  ///< Tracks total number of "test case" evaluations (across all organisms since beginning of run)

  void Setup();
  void SetupDiagnostic();
  void SetupSelection();
  void SetupEvaluation();
  void SetupMutator();
  void SetupDataCollection();

  template<typename DIAG_PROB>
  void SetupDiagnosticHelper();

  void SetupEvaluation_Cohort();
  void SetupEvaluation_DownSample();
  void SetupEvaluation_Full();

  void SetupSelection_Lexicase();
  void SetupSelection_Tournament();
  void SetupSelection_Truncation();
  void SetupSelection_None();
  void SetupSelection_Random();

  void InitializePopulation();

  void DoEvaluation();
  void DoSelection();
  void DoUpdate();

  void DoConfigSnapshot();
  void DoPopSnapshot();

public:

  DiagnosticsWorld(
    const DiagnosticsConfig& in_config
  ) :
    base_t("DiagnosticsWorld", false),
    config(in_config)
  {
    NewRandom(config.SEED()); // Set the seed in the world's random number generator.
    Setup();
  }

  ~DiagnosticsWorld() {
    // TODO - delete pointers!
    if (base_diagnostic != nullptr) base_diagnostic.Delete();
    if (selector != nullptr) selector.Delete();
  }

  void RunStep();
  void Run();

};

void DiagnosticsWorld::RunStep() {
  DoEvaluation();
  DoSelection();
  DoUpdate();
}

void DiagnosticsWorld::Run() {
  for (size_t u = 0; u <= config.MAX_GENS(); ++u) {
    RunStep();
  }
}

void DiagnosticsWorld::DoEvaluation() {
  emp_assert(org_aggregate_scores.size() == config.POP_SIZE());
  emp_assert(fit_fun_set.size() == config.POP_SIZE());
  emp_assert(org_test_scores.size() == config.POP_SIZE());
  emp_assert(org_test_evaluations.size() == config.POP_SIZE());
  // Assign test groupings (if any)
  assign_test_groupings();
  // Assign organism groupings (if any)
  assign_org_groupings();
  // Evaluate each organism
  for (size_t org_id = 0; org_id < GetSize(); ++org_id) {
    // Evaluation signal actions will:
    // - Translate organism genomes
    // - Update test scores, update aggregate score
    do_org_evaluation_sig.Trigger(org_id);
  }
}

void DiagnosticsWorld::DoSelection() {
  // Run selection-specific routine
  run_selection_routine();
  emp_assert(selected_parent_ids.size() == config.POP_SIZE());
  // each selected parent id reproduces
  for (size_t id : selected_parent_ids) {
    DoBirth(GetGenomeAt(id), id);
  }
}

void DiagnosticsWorld::DoUpdate() {
  // TODO
  // (1) Compute any per-generation statistics?
}


void DiagnosticsWorld::Setup() {
  std::cout << "--- Setting up DiagnosticsWorld ---" << std::endl;

  // Reset the world
  Reset();
  total_test_evaluations = 0;

  // Configure world to set organism ID on placement
  OnPlacement(
    [this](size_t pos) {
      auto& org = GetOrg(pos);
      org.SetPopID(pos);
      emp_assert(org.GetPopID() == pos);
    }
  );

  // TODO - Configure OnOffspringReady
  // OnOffspringReady(
  //   [this](org_t& org, size_t parent_pos) {
  //     // TODO
  //     // Mutate (make sure happens after systematics)
  //   }
  // );

  // Setup diagnostic problem
  SetupDiagnostic();

  // Setup evaluation
  SetupEvaluation();

  // Configure selection
  SetupSelection();

  // Setup mutation function
  SetupMutator();

  // Setup population structure
  SetPopStruct_Mixed(true);

  // TODO - disable (automatic?) mutations
  // Initialize population
  InitializePopulation();
  // TODO - print population check if what expected
  // TODO - enable (automatic?) mutations
  SetAutoMutate();

}

void DiagnosticsWorld::SetupDiagnostic() {
  std::cout << "Configuring diagnostic: " << config.DIAGNOSTIC() << std::endl;

  if (config.DIAGNOSTIC() == "exploitation-rate") {
    // SetupDiagnostic_ExploitationRate();
    SetupDiagnosticHelper<ExploitationRateDiagnostic>();
  } else if (config.DIAGNOSTIC() == "ordered-exploitation") {
    // SetupDiagnostic_OrderedExploitation();
    SetupDiagnosticHelper<OrderedExploitationDiagnostic>();
  } else if (config.DIAGNOSTIC() == "contradictory-objectives") {
    // SetupDiagnostic_ContradictoryObjectives();
    SetupDiagnosticHelper<ContradictoryObjectivesDiagnostic>();
  } else if (config.DIAGNOSTIC() == "multipath-exploration") {
    // SetupDiagnostic_MultipathExploration();
    SetupDiagnosticHelper<MultiPathExplorationDiagnostic>();
  } else {
    std::cout << "ERROR: UNKNOWN DIAGNOSTIC, " << config.DIAGNOSTIC() << std::endl;
    exit(-1);
  }
}

template<typename DIAG_PROB>
void DiagnosticsWorld::SetupDiagnosticHelper() {
  base_diagnostic = emp::NewPtr<DIAG_PROB>(
    config.TARGET(),
    config.CREDIT()
  );

  // Configure genome translation
  if (config.VALLEY_CROSSING()) {
    translate_genome_fun = [this](const genome_t& genome, phenotype_t& phen) {
      auto& diag = *(base_diagnostic.Cast<DIAG_PROB>());
      // Translate according to base diagnostic
      diag.Translate(genome, phen);
      // Apply valley crossing over phenotype
      phenotype_t valley_phen(valley_diagnostic.Translate(phen));
      emp_assert(valley_phen.size() == phen.size());
      std::copy(
        valley_phen.begin(),
        valley_phen.end(),
        phen.begin()
      );
    };
  } else {
    translate_genome_fun = [this](const genome_t& genome, phenotype_t& phen) {
      auto& diag = *(base_diagnostic.Cast<DIAG_PROB>());
      diag.Translate(genome, phen);
    };
  }

}

void DiagnosticsWorld::SetupMutator() {
  std::cout << "Configuring mutator" << std::endl;
  // TODO - Need to set things up to automutate?
  SetMutFun(
    [this](org_t& org, emp::Random& random) {
      // number of mutations and solution genome
      size_t mcnt = 0;
      genome_t& genome = org.GetGenome();

      // quick checks
      emp_assert(genome.size() == config.DIAGNOSTIC_DIMENSIONALITY());
      emp_assert(config.TARGET() > 0);

      for (size_t i = 0; i < genome.size(); ++i) {
        // if we do a mutation at this objective
        if (random.P(config.MUTATE_PER_SITE_RATE())) {
          const double mut = random.GetRandNormal(config.MUTATE_MEAN(), config.MUTATE_STD());
          if ( config.TARGET() < (genome[i] + mut) ) {
            // Rebound
            genome[i] = config.TARGET() - (genome[i] + mut - config.TARGET());
          } else if ( (genome[i] + mut) < config.GENE_LOWER_BND() ) {
            // Rebound
            genome[i] = std::abs(genome[i] + mut) + config.GENE_LOWER_BND();
          } else {
            // Add mutation
            genome[i] = genome[i] + mut;
          }
          ++mcnt;
        }
      }
      return mcnt;
    }
  );
}

void DiagnosticsWorld::SetupEvaluation() {
  std::cout << "Configuring evaluation (mode: " << config.EVAL_MODE() << ")" << std::endl;
  // Total tests is equal to diagnostic dimensionality.
  total_tests = config.DIAGNOSTIC_DIMENSIONALITY();
  // Allocate space for tracking organism aggregate scores
  org_aggregate_scores.clear();
  org_aggregate_scores.resize(config.POP_SIZE(), 0.0);
  // Allocate space for tracking organism test scores
  org_test_scores.clear();
  org_test_scores.resize(
    config.POP_SIZE(),
    emp::vector<double>(total_tests, 0.0)
  );
  // Allocate space for tracking organism test evaluations
  org_test_evaluations.clear();
  org_test_evaluations.resize(
    config.POP_SIZE(),
    emp::vector<bool>(total_tests, false)
  );
  // Initialize all possible test ids
  possible_test_ids.resize(total_tests, 0);
  std::iota(
    possible_test_ids.begin(),
    possible_test_ids.end(),
    0
  );
  // Initialize all possible population ids
  possible_pop_ids.resize(config.POP_SIZE(), 0);
  std::iota(
    possible_pop_ids.begin(),
    possible_pop_ids.end(),
    0
  );
  // Initialize org group ids
  org_group_ids.clear();
  org_group_ids.resize(config.POP_SIZE(), 0);

  // Clear all actions associated with organism evaluation.
  do_org_evaluation_sig.Clear();
  // First: translate the organism's genome, compute trait info
  do_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      emp_assert(org_id < this->GetSize());
      auto& org = this->GetOrg(org_id);
      org.TranslateGenome(translate_genome_fun);
      org.CalcOptimalTraits(config.TARGET(), config.ACCURACY());
    }
  );
  // Next: reset world organism info
  do_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      // Set all evaluated to false
      org_aggregate_scores[org_id] = 0.0;
      std::fill(
        org_test_scores[org_id].begin(),
        org_test_scores[org_id].end(),
        0.0
      );
      std::fill(
        org_test_evaluations[org_id].begin(),
        org_test_evaluations[org_id].end(),
        false
      );
    }
  );

  // Setup the default estimate_test_score functionality
  // TODO - setup configurable fitness estimation
  estimate_test_score = [this](size_t org_id, size_t test_id) {
    return org_test_scores[org_id][test_id];
  };

  // Configure the aggregate score functions
  agg_score_fun_set.clear();
  for (size_t org_id = 0; org_id < config.POP_SIZE(); ++org_id) {
    agg_score_fun_set.emplace_back(
      [this, org_id]() {
        emp_assert(org_id < org_aggregate_scores.size());
        return org_aggregate_scores[org_id];
      }
    );
  }

  // Configure the fitness functions (per-organism, per-test)
  fit_fun_set.clear();
  fit_fun_set.resize(config.POP_SIZE(), emp::vector<std::function<double(void)>>(0));
  for (size_t org_id = 0; org_id < config.POP_SIZE(); ++org_id) {
    for (size_t test_id = 0; test_id < total_tests; ++test_id) {
      fit_fun_set[org_id].emplace_back(
        [this, org_id, test_id]() {
          emp_assert(org_id < org_test_evaluations.size());
          emp_assert(org_id < org_test_scores.size());
          emp_assert(test_id < org_test_evaluations[org_id].size());
          emp_assert(test_id < org_test_scores[org_id].size());
          if (org_test_evaluations[org_id][test_id]) {
            return org_test_scores[org_id][test_id];
          } else {
            return estimate_test_score(org_id, test_id);
          }
        }
      );
    }

  }

  // OLD way of configuring fit_fun_set
  // for (size_t test_id = 0; test_id < total_tests; ++test_id) {
  //   fit_fun_set.emplace_back(
  //     [this, test_id](const org_t& org) {
  //       const size_t org_id = org.GetPopID();
  //       emp_assert(org_id < org_test_evaluations.size());
  //       emp_assert(org_id < org_test_scores.size());
  //       emp_assert(test_id < org_test_evaluations[org_id].size());
  //       emp_assert(test_id < org_test_scores[org_id].size());
  //       if (org_test_evaluations[org_id][test_id]) {
  //         return org_test_scores[org_id][test_id];
  //       } else {
  //         return estimate_test_score(org_id, test_id);
  //       }
  //     }
  //   );
  // }

  // full vs cohort vs down-sample
  if (config.EVAL_MODE() == "full") {
    SetupEvaluation_Full();
  } else if (config.EVAL_MODE() == "cohort") {
    SetupEvaluation_Cohort();
  } else if (config.EVAL_MODE() == "down-sample") {
    SetupEvaluation_DownSample();
  } else {
    std::cout << "Unknown EVAL_MODE: " << config.EVAL_MODE() << std::endl;
    exit(-1);
  }

}

void DiagnosticsWorld::SetupEvaluation_Cohort() {
  std::cout << "Configuring evaluation mode: cohort" << std::endl;
  emp_assert(config.NUM_COHORTS() > 0);
  emp_assert(total_tests > 0);
  emp_assert(config.POP_SIZE() > 0);
  const size_t num_cohorts = config.NUM_COHORTS();
  // Compute test cohort sizes
  const size_t base_test_cohort_size = (size_t)(total_tests / num_cohorts);
  size_t leftover_tests = total_tests - (base_test_cohort_size * num_cohorts);
  // Compute organism cohort sizes
  const size_t base_org_cohort_size = (size_t)(config.POP_SIZE() / num_cohorts);
  size_t leftover_orgs = config.POP_SIZE() - (base_org_cohort_size * num_cohorts);

  std::cout << "num_cohorts = " << num_cohorts << std::endl;
  std::cout << "base_test_cohort_size = " << base_test_cohort_size << std::endl;
  std::cout << "leftover_tests = " << leftover_tests << std::endl;

  // Initialize test groupings
  test_groupings.resize(num_cohorts);
  for (size_t i = 0; i < test_groupings.size(); ++i) {
    size_t group_size = base_test_cohort_size;
    if (leftover_tests > 0) {
      ++group_size;
      --leftover_tests;
    }
    auto& test_group = test_groupings[i];
    test_group.group_id = i;
    test_group.Resize(group_size, 0);
    std::cout << "  Test group " << i << " size: " << test_group.member_ids.size() << std::endl;
  }
  emp_assert(leftover_tests == 0);

  std::cout << "base_org_cohort_size = " << base_org_cohort_size << std::endl;
  std::cout << "leftover_orgs = " << leftover_orgs << std::endl;

  // Initialize org groupings
  org_groupings.resize(num_cohorts);
  for (size_t i = 0; i < org_groupings.size(); ++i) {
    size_t group_size = base_org_cohort_size;
    if (leftover_orgs > 0) {
      ++group_size;
      --leftover_orgs;
    }
    auto& org_group = org_groupings[i];
    org_group.group_id = i;
    org_group.Resize(group_size, 0);
    std::cout << "  Org group " << i << " size: " << org_group.member_ids.size() << std::endl;
  }
  emp_assert(leftover_orgs == 0);


  // Setup function to assign test/org groupings
  assign_test_groupings = [this]() {
    // Shuffle all possible test ids
    emp::Shuffle(*random_ptr, possible_test_ids);
    // Assign to cohorts in shuffled order
    size_t cur_pos = 0;
    for (size_t cohort_id = 0; cohort_id < test_groupings.size(); ++cohort_id) {
      auto& cohort = test_groupings[cohort_id];
      emp_assert(cohort.member_ids.size() == cohort.group_size);
      for (size_t test_i = 0; test_i < cohort.group_size; ++test_i) {
        cohort.member_ids[test_i] = possible_test_ids[cur_pos];
        ++cur_pos;
      }
    }
    emp_assert(cur_pos == total_tests);
  };

  // Setup function to assign orgnism groupings
  assign_org_groupings = [this]() {
    // Shuffle all possible test ids
    emp::Shuffle(*random_ptr, possible_pop_ids);
    // Assign to cohorts in shuffled order
    size_t cur_pos = 0;
    for (size_t cohort_id = 0; cohort_id < org_groupings.size(); ++cohort_id) {
      auto& cohort = org_groupings[cohort_id];
      emp_assert(cohort.member_ids.size() == cohort.group_size);
      for (size_t member_i = 0; member_i < cohort.group_size; ++member_i) {
        const size_t pop_id = possible_pop_ids[cur_pos];
        emp_assert(GetOrg(pop_id).GetPopID() == pop_id);
        cohort.member_ids[member_i] = pop_id;
        org_group_ids[cur_pos] = cohort.group_id;
        ++cur_pos;
      }
    }
    emp_assert(cur_pos == config.POP_SIZE());
  };

  // Configure organism evaluation (in cohort context)
  do_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      emp_assert(org_id < this->GetSize());
      auto& org = this->GetOrg(org_id);
      const size_t group_id = org_group_ids[org_id];
      // Evaluate organism on all tests in appropriate group
      emp_assert(group_id < test_groupings.size());
      auto& test_grouping = test_groupings[group_id];
      auto& cohort_test_ids = test_grouping.member_ids;
      double aggregate_score = 0.0;
      for (size_t i = 0; i < cohort_test_ids.size(); ++i) {
        const size_t test_id = cohort_test_ids[i];
        emp_assert(org.IsEvaluated());
        emp_assert(test_id < org.GetPhenotype().size());
        // Update test score
        org_test_scores[org_id][test_id] = org.GetPhenotype()[test_id];
        aggregate_score += org.GetPhenotype()[test_id];
        // Update evaluated
        org_test_evaluations[org_id][test_id] = true;
        ++total_test_evaluations;
      }
      // Update aggregate score
      org_aggregate_scores[org_id] = aggregate_score;
    }
  );

  // TODO print out assigned groupings, make sure makes sense


}

void DiagnosticsWorld::SetupEvaluation_DownSample() {
  std::cout << "Configuring Evaluation mode: down-sample" << std::endl;
  emp_assert(config.TEST_DOWNSAMPLE_RATE() > 0);
  emp_assert(config.TEST_DOWNSAMPLE_RATE() <= 1.0);
  emp_assert(total_tests > 0);

  size_t sample_size = (size_t)(config.TEST_DOWNSAMPLE_RATE() * (double)total_tests);
  sample_size = (sample_size == 0) ? sample_size + 1 : sample_size;
  emp_assert(sample_size > 0);
  emp_assert(sample_size <= total_tests);

  std::cout << "sample_size = " << sample_size << std::endl;

  // Initialize test groupings to one group that will contain random sample
  test_groupings.resize(1);
  auto& test_group = test_groupings.back();
  test_group.group_id = 0;
  test_group.Resize(sample_size, 0);
  emp_assert(test_group.member_ids.size() == sample_size);

  // Setup function to assign random test cases to test group
  assign_test_groupings = [this, sample_size]() {
    emp_assert(test_groupings.size() == 1);
    // Suffle all possible test ids
    emp::Shuffle(*random_ptr, possible_test_ids);
    // NOTE - if wanted to allow over sampling (> num tests), could modify this code to support
    auto& test_group = test_groupings.back();
    emp_assert(test_group.member_ids.size() == sample_size);
    for (size_t i = 0; i < sample_size; ++i) {
      const size_t test_id = possible_test_ids[i];
      test_group.member_ids[i] = test_id;
    }
  };

  // Initialize org groupings to one group that will contain entire population
  org_groupings.resize(1);
  auto& org_group = org_groupings.back();
  org_group.group_id = 0;
  org_group.Resize(config.POP_SIZE(), 0);
  // Go ahead and initialize organism grouping to contain all organism ids
  std::iota(
    org_group.member_ids.begin(),
    org_group.member_ids.end(),
    0
  );
  // Setup function to assign organism groupings (should do nothing, no need to modify current grouping)
  assign_org_groupings = [this]() {
    emp_assert(org_groupings.size() == 1);
    emp_assert(org_groupings.back().member_ids.size() == config.POP_SIZE());
    /* Do nothing */
  };

  // Configure organism evaluation
  do_org_evaluation_sig.AddAction(
    [this, sample_size](size_t org_id) {
      emp_assert(org_id < GetSize());
      emp_assert(test_groupings.size() == 1);
      auto& org = GetOrg(org_id);
      auto& test_group = test_groupings.back();
      emp_assert(org.IsEvaluated());
      emp_assert(test_group.member_ids.size() == sample_size);
      double aggregate_score = 0.0;
      for (size_t i = 0; i < sample_size; ++i) {
        const size_t test_id = test_group.member_ids[i];
        emp_assert(test_id < org.GetPhenotype().size());
        const double test_score = org.GetPhenotype()[test_id];
        // Update test score, aggregate score
        org_test_scores[org_id][test_id] = test_score;
        aggregate_score += test_score;
        // Update evaluated
        org_test_evaluations[org_id][test_id] = true;
        ++total_test_evaluations;
      }
      // Update aggregate score
      org_aggregate_scores[org_id] = aggregate_score;
    }
  );

  // TODO - print out assigned down sample, test to make sure it works

}

void DiagnosticsWorld::SetupEvaluation_Full() {
  std::cout << "Configuring evaluation mode: full" << std::endl;
  emp_assert(total_tests > 0);

  // Initialize the test groupings with one group that holds all tests
  test_groupings.resize(1);
  auto& test_group = test_groupings.back();
  test_group.group_id = 0;
  test_group.Resize(total_tests, 0);
  std::iota(
    test_group.member_ids.begin(),
    test_group.member_ids.end(),
    0
  );
  // Setup function to assign test cases to test group (should do nothing to change test group)
  assign_test_groupings = [this]() {
    emp_assert(test_groupings.size() == 1);
    emp_assert(test_groupings.back().member_ids.size() == total_tests);
    /* Do nothing */
  };

  // Initialize org groupings to one group containing all organisms
  org_groupings.resize(1);
  auto& org_group = org_groupings.back();
  org_group.group_id = 0;
  org_group.Resize(config.POP_SIZE(), 0);
  // Go ahead and initialize organism grouping to contain all organism ids
  std::iota(
    org_group.member_ids.begin(),
    org_group.member_ids.end(),
    0
  );
  // Setup function to assign organism groupings (should do nothing, no need to modify current grouping)
  assign_org_groupings = [this]() {
    emp_assert(org_groupings.size() == 1);
    emp_assert(org_groupings.back().member_ids.size() == config.POP_SIZE());
    /* Do nothing */
  };

  // Configure organism evaluation
  do_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      emp_assert(org_id < GetSize());
      emp_assert(test_groupings.size() == 1);
      auto& org = GetOrg(org_id);
      auto& test_group = test_groupings.back();
      emp_assert(org.IsEvaluated());
      emp_assert(test_group.member_ids.size() == total_tests);
      double aggregate_score = 0.0;
      for (size_t test_id = 0; test_id < total_tests; ++test_id) {
        emp_assert(test_id < org.GetPhenotype().size());
        const double test_score = org.GetPhenotype()[test_id];
        // Update test score, aggregate score
        org_test_scores[org_id][test_id] = test_score;
        aggregate_score += test_score;
        // Update evaluated
        org_test_evaluations[org_id][test_id] = true;
        ++total_test_evaluations;
      }
      // Update aggregate score
      org_aggregate_scores[org_id] = aggregate_score;
    }
  );

  // TODO - test full evaluation to make sure it works

}

void DiagnosticsWorld::SetupSelection() {
  std::cout << "Configuring parent selection routine" << std::endl;
  emp_assert(selector == nullptr);

  // Configure selection routine
  run_selection_routine = [this]() {
    // Resize parent ids to hold pop_size parents
    selected_parent_ids.resize(config.POP_SIZE(), 0);
    emp_assert(test_groupings.size() == org_groupings.size());
    const size_t num_groups = org_groupings.size();
    // For each grouping, select a number of parents equal to group size
    size_t num_selected = 0;
    for (size_t group_id = 0; group_id < num_groups; ++group_id) {
      auto& org_group = org_groupings[group_id];
      auto& test_group = test_groupings[group_id];
      const size_t n = org_group.GetSize();
      auto& selected = selection_fun(
        n,
        org_group.member_ids,
        test_group.member_ids
      );
      emp_assert(selected.size() == n);
      emp_assert(n + num_selected <= selected_parent_ids.size());
      std::copy(
        selected.begin(),
        selected.end(),
        selected_parent_ids.begin() + num_selected // TODO - check if this actually works!
      );
      num_selected += n;
    }
    // TODO - check that sets of selected ids correctly stored in selected_parent_ids
  };

  if (config.SELECTION() == "lexicase" ) {
    SetupSelection_Lexicase();
  } else if (config.SELECTION() == "tournament" ) {
    SetupSelection_Tournament();
  } else if (config.SELECTION() == "truncation" ) {
    SetupSelection_Truncation();
  } else if (config.SELECTION() == "none" ) {
    SetupSelection_None();
  } else if (config.SELECTION() == "random" ) {
    SetupSelection_Random();
  } else {
    std::cout << "Unknown selection scheme: " << config.SELECTION() << std::endl;
    exit(-1);
  }


}

void DiagnosticsWorld::SetupSelection_Lexicase() {
  // TODO - test lexicase
  selector = emp::NewPtr<selection::LexicaseSelect>(
    fit_fun_set,
    *random_ptr
  );

  selection_fun = [this](
    size_t n,
    const emp::vector<size_t>& org_group,
    const emp::vector<size_t>& test_group
  ) -> emp::vector<size_t>& {
    // Cast selector to lexicase selection
    auto& sel = *(selector.Cast<selection::LexicaseSelect>());
    return sel(n, org_group, test_group);
  };

}

void DiagnosticsWorld::SetupSelection_Tournament() {
  selector = emp::NewPtr<selection::TournamentSelect>(
    agg_score_fun_set,
    *random_ptr,
    config.TOURNAMENT_SIZE()
  );

  selection_fun = [this](
    size_t n,
    const emp::vector<size_t>& org_group,
    const emp::vector<size_t>& test_group
  ) -> emp::vector<size_t>& {
    // Cast selector to lexicase selection
    auto& sel = *(selector.Cast<selection::TournamentSelect>());
    return sel(n, org_group);
  };
}

void DiagnosticsWorld::SetupSelection_Truncation() {
  // TODO - setup truncation
  emp_assert(false);
}

void DiagnosticsWorld::SetupSelection_None() {
  // TODO - setup lexicase
  emp_assert(false);
}

void DiagnosticsWorld::SetupSelection_Random() {
  // TODO - setup lexicase
  emp_assert(false);
}


void DiagnosticsWorld::InitializePopulation() {

  std::cout << "Initializing population" << std::endl;

  if (config.INIT_POP_RAND()) {
    for (size_t i = 0; i < config.POP_SIZE(); ++i) {
      Inject(
        emp::RandomDoubleVector(
          *random_ptr,
          config.DIAGNOSTIC_DIMENSIONALITY(),
          config.GENE_LOWER_BND(),
          config.GENE_UPPER_BND()
        ),
        1
      );
    }
  } else {
    org_t default_org(config.DIAGNOSTIC_DIMENSIONALITY(), 0.0);
    Inject(default_org.GetGenome(), config.POP_SIZE());
  }

}


}