/**
 * Code adapted from https://github.com/jgh9094/ECJ-2022-suite-of-diagnostics-for-selection-schemes/ by Jose Hernandez
*/
#pragma once

#include <iostream>

#include "emp/Evolve/World.hpp"

#include "phylogeny/Phylogeny.hpp"

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

  // using eval_org_fun_t = std::function<void(org_t&)>;
  // using translate_genome_fun_t = std::function<void(const genome_t&, phenotype_t&)>;

  // struct TestGroup {
  //   emp::vector<size_t> test_ids;
  // };

  // struct PopGroup {
  //   emp::vector<size_t> pop_ids;
  // };

  struct Grouping {
    size_t group_id = 0;
    size_t group_size = 0;
    emp::vector<size_t> member_ids;

    void Resize(size_t size, size_t value=0) {
      group_size = size;
      member_ids.resize(group_size, value);
    }
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
  emp::vector< std::function<double(const org_t&)> > fit_fun_set; ///< One function for every possible test case.

  emp::vector< emp::vector<double> > org_test_scores; ///< Test scores for each organism
  emp::vector< emp::vector<bool> > org_test_evaluations; ///< Which test cases has each organism been evaluated on?

  emp::vector<size_t> possible_test_ids;
  emp::vector<size_t> possible_pop_ids;

  emp::vector<Grouping> test_groupings; ///< Groupings of tests
  emp::vector<Grouping> org_groupings;  ///< Groupings of organisms
  emp::vector<size_t> org_group_ids;    ///< Group ID for each organism

  std::function<void()> assign_test_groupings;
  std::function<void()> assign_org_groupings;

  std::function<void()> do_sample_tests;

  // emp::Ptr<BaseSelect> selector;
  emp::vector<size_t> selected_parent_ids;


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
  }

  void RunStep();
  void Run();

};

void DiagnosticsWorld::Setup() {
  std::cout << "--- Setting up DiagnosticsWorld ---" << std::endl;

  // Reset the world
  Reset();

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

  // // Configure organism evaluation
  // evaluate_org_fun = [this](org_t& org) {
  //   // Translate organism genome
  //   org.TranslateGenome(translate_genome_fun);
  //   // Calculate optimal traits
  //   org.CalcOptimalTraits(config.TARGET(), config.ACCURACY());
  // };
}

void DiagnosticsWorld::SetupMutator() {
  std::cout << "Configuring mutator" << std::endl;
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
          } else if ( (genome[i] + mut) < config.LOWER_BND() ) {
            // Rebound
            genome[i] = std::abs(genome[i] + mut) + config.LOWER_BND();
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

  // TODO - fit_fun_set

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

  // fit_fun_set
  // test_groupings
  // org_groupings

}

void DiagnosticsWorld::SetupEvaluation_Cohort() {
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
        cohort.member_ids[member_i] = possible_pop_ids[cur_pos];
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
      }
      // Update aggregate score
      org_aggregate_scores[org_id] = aggregate_score;
    }
  );

  // TODO print out assigned groupings, make sure makes sense


}

void DiagnosticsWorld::SetupEvaluation_DownSample() {
  // TODO
}

void DiagnosticsWorld::SetupEvaluation_Full() {
  // TODO
}

void DiagnosticsWorld::SetupSelection() {
  std::cout << "Configuring parent selection routine" << std::endl;

}


}