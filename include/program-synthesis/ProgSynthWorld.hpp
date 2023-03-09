#pragma once

#include <string>
#include <utility>
#include <algorithm>
#include <functional>

#include "emp/Evolve/World.hpp"
#include "emp/Evolve/Systematics.hpp"
#include "emp/bits/BitSet.hpp"
#include "emp/matching/MatchBin.hpp"
#include "emp/base/Ptr.hpp"
#include "emp/control/Signal.hpp"

// #include "sgp/cpu/linprg/LinearProgram.hpp"
// #include "sgp/cpu/LinearProgramCPU.hpp"
#include "sgp/cpu/lfunprg/LinearFunctionsProgram.hpp"
#include "sgp/cpu/LinearFunctionsProgramCPU.hpp"
#include "sgp/cpu/mem/BasicMemoryModel.hpp"
#include "sgp/inst/InstructionLibrary.hpp"
#include "sgp/EventLibrary.hpp"
#include "sgp/inst/lfpbm/InstructionAdder.hpp"

#include "../phylogeny/phylogeny_utils.hpp"
#include "../utility/Grouping.hpp"
#include "../selection/SelectionSchemes.hpp"

#include "ProgSynthConfig.hpp"
#include "ProgSynthOrg.hpp"
#include "Event.hpp"
#include "ProblemManager.hpp"
#include "ProgSynthHardware.hpp"
#include "MutatorLinearFunctionsProgram.hpp"

// TODO - re-organize problem manager <==> world interactions to use world signals
// i.e., pass the world to the problem manager configure, allow it to wire up functions to OnXSetup signals.

// TODO - alternatively, give the problem manager a reference to the world and
//        implement any necessary accessors

namespace psynth {

// TODO - move these into internal namespace
constexpr size_t TAG_SIZE = 32;
constexpr size_t FUNC_NUM_TAGS = 1;
constexpr size_t INST_TAG_CNT = 1;
constexpr size_t INST_ARG_CNT = 3;
using tag_t = emp::BitSet<TAG_SIZE>;
using inst_arg_t = int;
using program_t = sgp::cpu::lfunprg::LinearFunctionsProgram<tag_t, inst_arg_t>;
using ORGANISM_T = ProgSynthOrg<program_t>;
using MEMORY_MODEL_T = sgp::cpu::mem::BasicMemoryModel;
using MATCHBIN_T = emp::MatchBin<
  size_t,
  emp::HammingMetric<TAG_SIZE>,
  emp::RankedSelector<>,
  emp::NopRegulator
>;

class ProgSynthWorld : public emp::World<ORGANISM_T> {
public:
  // --- Type aliases --
  using org_t = ORGANISM_T;
  using base_t = emp::World<org_t>;
  using genome_t = typename org_t::genome_t;
  using phenotype_t = typename org_t::phenotype_t;
  using hw_memory_model_t = MEMORY_MODEL_T;
  using hw_matchbin_t = MATCHBIN_T;
  using inst_t = typename program_t::inst_t;
  using hardware_t = sgp::cpu::LinearFunctionsProgramCPU<
    hw_memory_model_t,
    inst_arg_t,
    hw_matchbin_t,
    ProgSynthHardwareComponent<tag_t>
  >;
  using inst_lib_t = sgp::inst::InstructionLibrary<hardware_t, inst_t>;
  using event_lib_t = sgp::EventLibrary<hardware_t>;
  using base_event_t = typename event_lib_t::event_t;
  using mutator_t = MutatorLinearFunctionsProgram<hardware_t, tag_t, inst_arg_t>;
  using selection_fun_t = std::function<
    emp::vector<size_t>&(
      size_t,
      const emp::vector<size_t>&,
      const emp::vector<size_t>&
    )
  >;
  using systematics_t = emp::Systematics<
    org_t,
    genome_t,
    phylo::phenotype_info
  >;
  using taxon_info_t = phylo::phenotype_info;
  using taxon_t = typename systematics_t::taxon_t;

  using config_t = ProgSynthConfig;

  struct SelectedStatistics {
    // TODO
  };

protected:
  const config_t& config;

  emp::Ptr<hardware_t> eval_hardware = nullptr;
  inst_lib_t inst_lib;
  event_lib_t event_lib;
  emp::Ptr<mutator_t> mutator = nullptr;

  ProblemManager<hardware_t> problem_manager;
  // size_t event_id_input_sig = 0;
  size_t event_id_numeric_input_sig = 0;

  size_t total_tests = 0;
  size_t total_test_evaluations = 0; ///< Tracks the total number of "test case" evaluations across all organisms since the beginning of the run.
  bool found_solution = false;

  std::function<bool(void)> stop_run;

  emp::vector<double> org_aggregate_scores;
  emp::vector<bool> pop_test_coverage;

  emp::Signal<void(size_t)> do_org_evaluation_sig;

  emp::Signal<void(org_t&, size_t, bool)> begin_program_test_sig;
  emp::Signal<void(org_t&, size_t)> do_program_test_sig;
  emp::Signal<void(org_t&, size_t)> end_program_test_sig;

  emp::vector<
    emp::vector< std::function<double(void)> >
  > fit_fun_set;      ///< Per-organism, per-test
  emp::vector<
    std::function<double(void)>
  > agg_score_fun_set; ///< Per-organism, aggregate score

  std::function<double(size_t, size_t)> estimate_test_score;

  emp::vector< emp::vector<double> > org_test_scores;   ///< Test scores for each organism
  emp::vector< emp::vector<bool> > org_test_evaluations; ///< Which test cases has each organism been evaluated on?

  utils::GroupManager org_groupings;
  utils::GroupManager test_groupings;

  emp::Ptr<selection::BaseSelect> selector;
  emp::vector<size_t> selected_parent_ids;

  std::function<void(void)> run_selection_routine;
  selection_fun_t selector_fun;

  emp::Ptr<systematics_t> systematics_ptr;

  std::string output_dir;

  // -- Data files --
  emp::Ptr<emp::DataFile> summary_file_ptr;
  emp::Ptr<emp::DataFile> phylodiversity_file_ptr;
  emp::Ptr<emp::DataFile> elite_file_ptr;

  void Setup();
  void SetupProblem();
  void SetupSelection();
  void SetupEvaluation();
  void SetupFitFunEstimator();
  void SetupMutator();
  void SetupPhylogenyTracking();
  void SetupDataCollection();
  void SetupStoppingCondition();
  void SetupVirtualHardware();
  void SetupInstructionLibrary();
  void SetupEventLibrary();

  void SetupEvaluation_Full();
  void SetupEvaluation_Cohort();
  void SetupEvaluation_DownSample();

  void InitializePopulation();
  void DoEvaluation();
  void DoSelection();
  void DoUpdate();

  void SnapshotConfig();

public:

  ProgSynthWorld(
    const config_t& in_config
  ) :
    base_t("ProgSynthWorld", false),
    config(in_config)
  {
    NewRandom(config.SEED());
    Setup();
  }

  ~ProgSynthWorld() {
    // TODO
    if (eval_hardware != nullptr) { eval_hardware.Delete(); }
    if (mutator != nullptr) { mutator.Delete(); }
  }

  void RunStep();
  void Run();

  const config_t& GetConfig() const { return config; }

};

void ProgSynthWorld::RunStep() {
  // TODO - DoEvaluation();
  // TODO - DoSelection();
  // TODO - DoUpdate();
}

void ProgSynthWorld::DoEvaluation() {
  // TODO
}

void ProgSynthWorld::Setup() {
  std::cout << "--- Setting up ProgSynth ---" << std::endl;

  // Reset the world
  Reset();
  total_test_evaluations = 0;
  found_solution = false;

  // Setup the population structure
  SetPopStruct_Mixed(true);


  // Configure world to set organism ID on placement
  OnPlacement(
    [this](size_t pos) {
      auto& org = GetOrg(pos);
      org.SetPopID(pos);
      org.GetPhenotype().Reset(total_tests);
      emp_assert(org.GetPopID() == pos);
    }
  );

  // Setup the program synthesis problem
  SetupProblem();

  // Setup the instruction library
  SetupInstructionLibrary();

  // Setup the event library
  SetupEventLibrary();

  // Setup the virtual hardware used to evaluate programs
  SetupVirtualHardware();

  // Setup the program mutator
  SetupMutator();

  // Setup evaluation
  SetupEvaluation();

  // Setup selection
  // SetupSelection();


  // TODO - SetAutoMutate!
}

void ProgSynthWorld::SetupProblem() {
  std::cout << "Setting up problem." << std::endl;
  if (!problem_manager.IsValidProblem(config.PROBLEM())) {
    std::cout << "Unknown problem: " << config.PROBLEM() << std::endl;
    exit(-1);
  }
  problem_manager.ConfigureProblem(config.PROBLEM());
  problem_manager.LoadTrainingSet(config.TRAINING_SET_PATH());
  problem_manager.LoadTestingSet(config.TESTING_SET_PATH());
  // TODO - be sure to output training / testing set sizes in cfg snapshot
  std::cout << "  - Loaded training set size: " << problem_manager.GetTrainingSetSize() << std::endl;
  std::cout << "  - Loaded testing set size: " << problem_manager.GetTestingSetSize() << std::endl;
}

void ProgSynthWorld::SetupInstructionLibrary() {
  std::cout << "Setting up instruction library." << std::endl;
  // Reset instruction library
  inst_lib.Clear();
  sgp::inst::lfpbm::InstructionAdder<hardware_t> inst_adder;
  // Add default instructions
  inst_adder.AddAllDefaultInstructions(inst_lib);
  // Add problem-specific instructions
  problem_manager.AddProblemInstructions(inst_lib);
  // TODO - snapshot instruction set
}

void ProgSynthWorld::SetupEventLibrary() {
  std::cout << "Setting up event library." << std::endl;
  event_lib.Clear();
  // Add default event set
  event_id_numeric_input_sig = event_lib.AddEvent(
    "NumericInputSignal",
    [this](hardware_t& hw, const base_event_t& e) {
      const NumericMessageEvent<TAG_SIZE>& event = static_cast<const NumericMessageEvent<TAG_SIZE>&>(e);
      auto thread_id = hw.SpawnThreadWithTag(event.GetTag());
      if (thread_id && event.GetData().size()) {
        // If message resulted in thread being spawned, load message into local working space.
        auto& thread = hw.GetThread(thread_id.value());
        // Wait, wait. Does this thread have calls on the call stack?
        if (thread.GetExecState().call_stack.size()) {
          auto& call_state = thread.GetExecState().GetTopCallState();
          auto& mem_state = call_state.GetMemory();
          for (auto mem : event.GetData()) {
            mem_state.SetWorking(mem.first, mem.second);
          }
        }
      }
    }
  );
  // Configure problem-specific events
  problem_manager.AddProblemEvents(event_lib);
}

void ProgSynthWorld::SetupVirtualHardware() {
  std::cout << "Setting up virtual hardware." << std::endl;
  // TODO - Implement ability to run hardware in parallel
  if (eval_hardware == nullptr) {
    eval_hardware = emp::NewPtr<hardware_t>(*random_ptr, inst_lib, event_lib);
  }
  // Configure the SGP CPU
  eval_hardware->Reset();
  eval_hardware->SetActiveThreadLimit(config.MAX_ACTIVE_THREAD_CNT());
  eval_hardware->SetThreadCapacity(config.MAX_THREAD_CAPACITY());
  // Configure input tag to all 0s
  tag_t input_tag;
  input_tag.Clear();
  eval_hardware->GetCustomComponent().SetInputTag(input_tag);
  // Configure problem-specific hardware component.
  problem_manager.AddProblemHardware(*eval_hardware);
  // Hardware should be in a valid thread state after configuration.
  emp_assert(eval_hardware->ValidateThreadState());
}

void ProgSynthWorld::SetupMutator() {
  std::cout << "Setting up program mutator." << std::endl;
  if (mutator != nullptr) {
    mutator.Delete();
  }
  mutator = emp::NewPtr<mutator_t>(inst_lib);
  mutator->ResetLastMutationTracker();
  // Set program constraints
  mutator->SetProgFunctionCntRange(
    {config.PRG_MIN_FUNC_CNT(), config.PRG_MAX_FUNC_CNT()}
  );
  mutator->SetProgFunctionInstCntRange(
    {config.PRG_MIN_FUNC_INST_CNT(), config.PRG_MAX_FUNC_INST_CNT()}
  );
  mutator->SetProgInstArgValueRange(
    {config.PRG_INST_MIN_ARG_VAL(), config.PRG_INST_MAX_ARG_VAL()}
  );
  const size_t total_inst_limit = 2 * config.PRG_MAX_FUNC_INST_CNT() * config.PRG_MAX_FUNC_CNT();
  mutator->SetTotalInstLimit(total_inst_limit);
  mutator->SetFuncNumTags(FUNC_NUM_TAGS);
  mutator->SetInstNumTags(INST_TAG_CNT);
  mutator->SetInstNumArgs(INST_ARG_CNT);
  // Set mutation rates
  mutator->SetRateInstArgSub(config.MUT_RATE_INST_ARG_SUB());
  mutator->SetRateInstTagBF(config.MUT_RATE_INST_TAG_BF());
  mutator->SetRateInstSub(config.MUT_RATE_INST_SUB());
  mutator->SetRateInstIns(config.MUT_RATE_INST_INS());
  mutator->SetRateInstDel(config.MUT_RATE_INST_DEL());
  mutator->SetRateSeqSlip(config.MUT_RATE_SEQ_SLIP());
  mutator->SetRateFuncDup(config.MUT_RATE_FUNC_DUP());
  mutator->SetRateFuncDel(config.MUT_RATE_FUNC_DEL());
  mutator->SetRateFuncTagBF(config.MUT_RATE_FUNC_TAG_BF());
  mutator->SetRateInstTagSingleBF(config.MUT_RATE_INST_TAG_SINGLE_BF());
  mutator->SetRateFuncTagSingleBF(config.MUT_RATE_FUNC_TAG_SINGLE_BF());
  mutator->SetRateInstTagSeqRand(config.MUT_RATE_INST_TAG_SEQ_RAND());
  mutator->SetRateFuncTagSeqRand(config.MUT_RATE_FUNC_TAG_SEQ_RAND());

  // Set world mutation function
  SetMutFun(
    [this](org_t& org, emp::Random& rnd) {
      mutator->ResetLastMutationTracker(); // Reset mutator's recorded mutations.
      const size_t mut_cnt = mutator->ApplyAll(
        rnd,
        org.GetGenome().GetProgram()
      );
      return mut_cnt;
    }
  );
}

void ProgSynthWorld::SetupEvaluation() {
  std::cout << "Configuring evaluation (mode: " << config.EVAL_MODE() << ")" << std::endl;
  // Total tests is equal to number of tests we loaded into the testing set.
  total_tests = problem_manager.GetTestingSetSize();
  // Allocate space for tracking organism aggregate scores
  org_aggregate_scores.clear();
  org_aggregate_scores.resize(config.POP_SIZE(), 0.0);
  // Allocate space for tracking population-wide test coverage
  pop_test_coverage.clear();
  pop_test_coverage.resize(total_tests, false);
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
  // Setup organism group manager.
  emp::vector<size_t> possible_ids(config.POP_SIZE(), 0);
  std::iota(
    possible_ids.begin(),
    possible_ids.end(),
    0
  );
  org_groupings.SetPossibleIDs(possible_ids);
  std::cout << org_groupings.GetPossibleIDs() << std::endl;

  // Setup test group manager
  possible_ids.resize(total_tests, 0);
  std::iota(
    possible_ids.begin(),
    possible_ids.end(),
    0
  );
  test_groupings.SetPossibleIDs(possible_ids);
  std::cout << test_groupings.GetPossibleIDs() << std::endl;

  // Clear out all actions associated with organism evaluation.
  do_org_evaluation_sig.Clear();
  // TODO - configure organism evaluation
  // (1) Load program into hardware, run on all tests
  // (2) Update world's score tracking vectors
  // (3)

  begin_program_test_sig.AddAction(
    [this](org_t& org, size_t test_id, bool training) {
      eval_hardware->ResetMatchBin();      // Reset the matchbin between tests
      eval_hardware->ResetHardwareState(); // Reset hardware execution state information (global memory, threads, etc)
      // TODO - anything else needs to happen before test eval?
      // TODO - load test input
      problem_manager.InitCase(
        *eval_hardware,
        org,
        test_id,
        training
      );
    }
  );

  // TODO - end_program_test_sig

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
  // TODO - have test scores in one place? (organism or here)
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

  // Setup evaluation mode
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

  // Record taxon info after evaluation, but before summing aggregate scores
  // TODO

}

void ProgSynthWorld::SetupEvaluation_Full() {
  std::cout << "Configuring evaluation mode: full" << std::endl;
  emp_assert(total_tests > 0);

  // Initialize the test groupings with one group that holds all tests.
  test_groupings.SetSingleGroupMode();
  // Initialize the organism groupings with one group that holds all organisms.
  org_groupings.SetSingleGroupMode();
  emp_assert(test_groupings.GetNumGroups() == org_groupings.GetNumGroups());

  // Configure organism evaluation
  // TODO
  do_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      emp_assert(org_id < GetSize());
      emp_assert(test_groupings.GetNumGroups() == 1);
      auto& org = GetOrg(org_id);
      // TODO - this is common across evaluation modes?
      //        if so, move up
      // Reset phenotype
      phenotype_t& phen = org.GetPhenotype();
      phen.Reset(total_tests);
      // Load program onto evaluation hardware unit
      eval_hardware->SetProgram(org.GetGenome().GetProgram());

      const auto& test_group = test_groupings.GetGroup(0);
      const auto& test_ids = test_group.GetMembers();
      for (size_t i = 0; i < test_ids.size(); ++i) {
        const size_t test_id = test_ids[i]; // Get test id from group.
        // TODO - Run program on test.
        begin_program_test_sig.Trigger(org, test_id, true);
        do_program_test_sig.Trigger(org, test_id);
        end_program_test_sig.Trigger(org, test_id);
        emp_assert(false, "Have not implemented evaluation signals yet.");
        // TODO - Update scores / test
        org_test_evaluations[org_id][test_id] = true;
        ++total_test_evaluations;
      }
    }
  );

}

void ProgSynthWorld::SetupEvaluation_Cohort() {
  // TODO
}

void ProgSynthWorld::SetupEvaluation_DownSample() {
  // TODO
}

}