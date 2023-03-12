#pragma once

#include <string>
#include <utility>
#include <algorithm>
#include <functional>
#include <iostream>
#include <sys/stat.h>

#include "emp/Evolve/World.hpp"
#include "emp/Evolve/Systematics.hpp"
#include "emp/bits/BitSet.hpp"
#include "emp/matching/MatchBin.hpp"
#include "emp/base/Ptr.hpp"
#include "emp/control/Signal.hpp"

#include "sgp/cpu/lfunprg/LinearFunctionsProgram.hpp"
#include "sgp/cpu/LinearFunctionsProgramCPU.hpp"
#include "sgp/cpu/mem/BasicMemoryModel.hpp"
#include "sgp/inst/InstructionLibrary.hpp"
#include "sgp/EventLibrary.hpp"
#include "sgp/inst/lfpbm/InstructionAdder.hpp"

#include "../phylogeny/phylogeny_utils.hpp"
#include "../utility/Grouping.hpp"
#include "../utility/printing.hpp"
#include "../selection/SelectionSchemes.hpp"

#include "ProgSynthConfig.hpp"
#include "ProgSynthOrg.hpp"
#include "Event.hpp"
#include "ProblemManager.hpp"
#include "ProgSynthHardware.hpp"
#include "MutatorLinearFunctionsProgram.hpp"
#include "SelectedStatistics.hpp"

// TODO - implement program json output / input

// TODO - re-organize problem manager <==> world interactions to use world signals
// i.e., pass the world to the problem manager configure, allow it to wire up functions to OnXSetup signals.

// TODO - alternatively, give the problem manager a reference to the world and
//        implement any necessary accessors

namespace psynth {

namespace world_defs {

constexpr size_t TAG_SIZE = 32;
constexpr size_t FUNC_NUM_TAGS = 1;
constexpr size_t INST_TAG_CNT = 1;
constexpr size_t INST_ARG_CNT = 3;
using TAG_T = emp::BitSet<TAG_SIZE>;
using INST_ARG_T = int;
using PROGRAM_T = sgp::cpu::lfunprg::LinearFunctionsProgram<TAG_T, INST_ARG_T>;
using ORGANISM_T = ProgSynthOrg<PROGRAM_T>;
using MEMORY_MODEL_T = sgp::cpu::mem::BasicMemoryModel;
using MATCHBIN_T = emp::MatchBin<
  size_t,
  emp::HammingMetric<TAG_SIZE>,
  emp::RankedSelector<>,
  emp::NopRegulator
>;

}

class ProgSynthWorld : public emp::World<world_defs::ORGANISM_T> {
public:
  // --- Type aliases --
  using org_t = world_defs::ORGANISM_T;
  using base_t = emp::World<org_t>;
  using genome_t = typename org_t::genome_t;
  using phenotype_t = typename org_t::phenotype_t;
  using hw_memory_model_t = world_defs::MEMORY_MODEL_T;
  using hw_matchbin_t = world_defs::MATCHBIN_T;
  using program_t = world_defs::PROGRAM_T;
  using inst_t = typename program_t::inst_t;
  using inst_arg_t = world_defs::INST_ARG_T;
  using tag_t = world_defs::TAG_T;
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

  static constexpr size_t TAG_SIZE = world_defs::TAG_SIZE;
  static constexpr size_t FUNC_NUM_TAGS = world_defs::FUNC_NUM_TAGS;
  static constexpr size_t INST_TAG_CNT = world_defs::INST_TAG_CNT;
  static constexpr size_t INST_ARG_CNT = world_defs::INST_ARG_CNT;

protected:
  const config_t& config;

  bool world_configured = false;

  emp::Ptr<hardware_t> eval_hardware = nullptr;
  inst_lib_t inst_lib;
  event_lib_t event_lib;
  emp::Ptr<mutator_t> mutator = nullptr;

  ProblemManager<hardware_t> problem_manager;
  // size_t event_id_input_sig = 0;
  size_t event_id_numeric_input_sig = 0;

  size_t total_training_cases = 0;
  size_t total_test_evaluations = 0; ///< Tracks the total number of "test case" evaluations across all organisms since the beginning of the run.
  bool found_solution = false;

  std::function<bool(void)> stop_run;
  std::function<bool(void)> is_final_update;
  std::function<bool(size_t)> check_org_solution;

  emp::Signal<void(size_t)> begin_org_evaluation_sig;
  emp::Signal<void(size_t)> do_org_evaluation_sig;  // Overall evaluate organism process.
  emp::Signal<void(size_t)> end_org_evaluation_sig;

  emp::Signal<void(org_t&)> begin_program_eval_sig; // Triggered at beginning of program evaluation (program loaded on hardware, phenotype reset).

  emp::Signal<void(org_t&, size_t, bool)> begin_program_test_sig;
  emp::Signal<void(org_t&, size_t)> do_program_test_sig;
  emp::Signal<void(org_t&, size_t, bool)> end_program_test_sig;

  emp::vector<
    emp::vector< std::function<double(void)> >
  > fit_fun_set;      ///< Per-organism, per-test
  emp::vector<
    std::function<double(void)>
  > agg_score_fun_set; ///< Per-organism, aggregate score

  std::function<double(size_t, size_t)> estimate_test_score;

  emp::vector<bool> pop_training_coverage;
  emp::vector<double> org_aggregate_scores;
  emp::vector<size_t> org_training_coverage;  ///< Organism coverage of training cases
  emp::vector<size_t> org_num_training_cases; ///< Number of training cases organism has been evaluated against
  emp::vector< emp::vector<double> > org_training_scores;   ///< Test scores for each organism
  emp::vector< emp::vector<bool> > org_training_evaluations; ///< Which test cases has each organism been evaluated on?

  utils::GroupManager org_groupings;
  utils::GroupManager test_groupings;

  emp::vector<size_t> all_training_case_ids;
  emp::Ptr<selection::BaseSelect> selector = nullptr;
  emp::vector<size_t> selected_parent_ids;

  std::function<void(void)> run_selection_routine;
  selection_fun_t selection_fun;

  emp::Ptr<systematics_t> systematics_ptr = nullptr;


  // -- Output --
  std::string output_dir;

  emp::Ptr<emp::DataFile> summary_file_ptr = nullptr;
  emp::Ptr<emp::DataFile> phylodiversity_file_ptr = nullptr;
  emp::Ptr<emp::DataFile> elite_file_ptr = nullptr;

  SelectedStatistics selection_stats;

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

  void SetupSelection_Lexicase();
  void SetupSelection_Tournament();
  void SetupSelection_Truncation();
  void SetupSelection_None();
  void SetupSelection_Random();

  void SetupFitFunEstimator_None();
  void SetupFitFunEstimator_Ancestor();
  void SetupFitFunEstimator_Relative();

  void SetupStoppingCondition_Generations();
  void SetupStoppingCondition_Evaluations();

  void SetupDataCollection_Phylodiversity();
  void SetupDataCollection_Summary();
  void SetupDataCollection_Elite();

  void InitializePopulation();
  void InitializePopulation_Load();
  void InitializePopulation_Random();

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
    if (eval_hardware != nullptr) { eval_hardware.Delete(); }
    if (mutator != nullptr) { mutator.Delete(); }
    if (selector != nullptr) { selector.Delete(); }
    if (phylodiversity_file_ptr != nullptr) { phylodiversity_file_ptr.Delete(); }
    if (summary_file_ptr != nullptr) { summary_file_ptr.Delete(); }
    if (elite_file_ptr != nullptr) { elite_file_ptr.Delete(); }
  }

  void RunStep();
  void Run();

  const config_t& GetConfig() const { return config; }

  size_t GetNumTrainingCases() const {
    emp_assert(world_configured);
    return total_training_cases;
  }

};

void ProgSynthWorld::RunStep() {
  DoEvaluation();
  DoSelection();
  DoUpdate();
}

void ProgSynthWorld::Run() {
  while (!stop_run()) {
    RunStep();
  }
}

void ProgSynthWorld::DoEvaluation() {
  emp_assert(pop_training_coverage.size() == total_training_cases);
  emp_assert(org_aggregate_scores.size() == config.POP_SIZE());
  emp_assert(fit_fun_set.size() == config.POP_SIZE());
  emp_assert(org_training_scores.size() == config.POP_SIZE());
  emp_assert(org_training_evaluations.size() == config.POP_SIZE());
  emp_assert(org_training_coverage.size() == config.POP_SIZE());
  emp_assert(org_num_training_cases.size() == config.POP_SIZE());

  // TODO - reset true max fitness organism

  // Update test and organism groupings
  org_groupings.UpdateGroupings();
  test_groupings.UpdateGroupings();
  // Reset things...
  std::fill(
    pop_training_coverage.begin(),
    pop_training_coverage.end(),
    false
  );

  // Evaluate each organism
  for (size_t org_id = 0; org_id < GetSize(); ++org_id) {
    do_org_evaluation_sig.Trigger(org_id);
  }

  // Estimate aggregate fitness for each organism
  // Need to do this separately to ensure all possible descendants of current taxa have been evaluated
  // Otherwise, estimation can crash
  for (size_t org_id = 0; org_id < GetSize(); ++org_id) {
    double est_agg_score = 0.0;
    for (size_t test_id = 0; test_id < total_training_cases; ++test_id) {
      est_agg_score += fit_fun_set[org_id][test_id]();
    }
    org_aggregate_scores[org_id] = est_agg_score;
  }

}

void ProgSynthWorld::DoSelection() {
  // Run configured selection routine
  run_selection_routine();
  emp_assert(selected_parent_ids.size() == config.POP_SIZE());
  // Each selected parent id reproduces
  for (size_t id : selected_parent_ids) {
    DoBirth(GetGenomeAt(id), id);
  }
}

void ProgSynthWorld::DoUpdate() {
  // (1) Compute any per-generation statistics / intervals?
  emp_assert(config.PRINT_INTERVAL() > 0);
  const size_t cur_update = GetUpdate();
  const bool final_update = is_final_update();
  const bool print_interval = !(cur_update % config.PRINT_INTERVAL()) || final_update;
  const bool summary_data_interval = !(cur_update % config.OUTPUT_SUMMARY_DATA_INTERVAL()) || final_update;
  const bool snapshot_interval = !(cur_update % config.SNAPSHOT_INTERVAL()) || final_update;

  // (2) File output
  if (summary_data_interval) {
    // Update selection statistics
    selection_stats.Calculate(
      selected_parent_ids,
      *this
    );
    // TODO - Update files
    // summary_file_ptr->Update();
    // elite_file_ptr->Update();
    // phylodiversity_file_ptr->Update();
  }

  if (snapshot_interval) {
    // TODO - snapshot phylogeny
  }

  // (3) Print status
  if (print_interval) {
    std::cout << "update: " << GetUpdate() << "; ";
    // TODO - add fitness
    // std::cout << "best score (" << true_max_fit_org_id << "): " << GetOrg(true_max_fit_org_id).GetAggregateScore();
    std::cout << std::endl;
  }

  // (4) Update!
  Update();
}

void ProgSynthWorld::Setup() {
  std::cout << "--- Setting up ProgSynth ---" << std::endl;
  world_configured = false;

  // Reset the world
  Reset();
  total_test_evaluations = 0;
  found_solution = false;

  // Configure output directory path, create directory
  output_dir = config.OUTPUT_DIR();
  mkdir(output_dir.c_str(), ACCESSPERMS);
  if(output_dir.back() != '/') {
      output_dir += '/';
  }

  // Setup the population structure
  SetPopStruct_Mixed(true);


  // Configure world to set organism ID on placement
  OnPlacement(
    [this](size_t pos) {
      auto& org = GetOrg(pos);
      org.SetPopID(pos);
      org.GetPhenotype().Reset(total_training_cases);
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
  SetupSelection();
  // Setup fitness function estimator
  SetupFitFunEstimator();
  // Setup stopping condition
  SetupStoppingCondition();
  // Setup phylogeny tracking
  SetupPhylogenyTracking();
  // Setup data collection
  SetupDataCollection();
  // Initialize population!
  InitializePopulation();
  // SetAutoMutate!
  SetAutoMutate();
  // Output a snapshot of the run configuration
  SnapshotConfig();
  world_configured = true;
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
  inst_adder.AddAllDefaultInstructions(
    inst_lib,
    {"Fork", "Terminate"}
  );
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
  total_training_cases = problem_manager.GetTrainingSetSize();
  // Allocate space for tracking organism training coverage
  org_training_coverage.clear();
  org_training_coverage.resize(config.POP_SIZE(), 0);
  // Allocate space for tracking number of training cases an organism is evaluated against
  org_num_training_cases.clear();
  org_num_training_cases.resize(config.POP_SIZE(), 0);
  // Allocate space for tracking organism aggregate scores
  org_aggregate_scores.clear();
  org_aggregate_scores.resize(config.POP_SIZE(), 0.0);
  // Allocate space for tracking population-wide test coverage
  pop_training_coverage.clear();
  pop_training_coverage.resize(total_training_cases, false);
  // Allocate space for tracking organism test scores
  org_training_scores.clear();
  org_training_scores.resize(
    config.POP_SIZE(),
    emp::vector<double>(total_training_cases, 0.0)
  );
  // Allocate space for tracking organism test evaluations
  org_training_evaluations.clear();
  org_training_evaluations.resize(
    config.POP_SIZE(),
    emp::vector<bool>(total_training_cases, false)
  );
  // Setup organism group manager.
  emp::vector<size_t> possible_org_ids(config.POP_SIZE(), 0);
  std::iota(
    possible_org_ids.begin(),
    possible_org_ids.end(),
    0
  );
  org_groupings.SetPossibleIDs(possible_org_ids);
  // std::cout << "Possible org ids: " << org_groupings.GetPossibleIDs() << std::endl;

  // Setup test group manager
  all_training_case_ids.resize(total_training_cases, 0);
  std::iota(
    all_training_case_ids.begin(),
    all_training_case_ids.end(),
    0
  );
  test_groupings.SetPossibleIDs(all_training_case_ids);
  // std::cout << "Possible training case ids: " << test_groupings.GetPossibleIDs() << std::endl;

  // Clear out all actions associated with organism evaluation.
  begin_org_evaluation_sig.Clear();
  do_org_evaluation_sig.Clear();
  end_org_evaluation_sig.Clear();

  begin_program_eval_sig.Clear();
  begin_program_test_sig.Clear();
  do_program_test_sig.Clear();
  end_program_test_sig.Clear();

  // TODO - configure organism evaluation
  begin_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      // 0-out aggregate score
      org_aggregate_scores[org_id] = 0.0;
      // 0-out training coverage
      org_training_coverage[org_id] = 0;
      // 0-out num training cases evaluated against
      org_num_training_cases[org_id] = 0;
      // 0 out organism's training score
      std::fill(
        org_training_scores[org_id].begin(),
        org_training_scores[org_id].end(),
        0.0
      );
      // 0 out organism's evaluations
      std::fill(
        org_training_evaluations[org_id].begin(),
        org_training_evaluations[org_id].end(),
        false
      );
    }
  );

  end_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      // TODO - check if candidate for testing against testing set
      // - Need to check if num_passes == number of tests evaluated on
      // - If so, add to vector of ids to be tested whether they are solutions
      // - Allow each evaluation mode to implement the actual testing
      // TODO - record taxon information
    }
  );

  begin_program_eval_sig.AddAction(
    [this](org_t& org) {
      // Reset phenotype
      phenotype_t& phen = org.GetPhenotype();
      phen.Reset(total_training_cases);
      // Load program onto evaluation hardware unit
      eval_hardware->SetProgram(org.GetGenome().GetProgram());
    }
  );

  begin_program_test_sig.AddAction(
    [this](org_t& org, size_t test_id, bool training) {
      eval_hardware->ResetMatchBin();      // Reset the matchbin between tests
      eval_hardware->ResetHardwareState(); // Reset hardware execution state information (global memory, threads, etc)
      eval_hardware->GetCustomComponent().Reset(); // Reset custom component
      // TODO - anything else needs to happen before test eval?
      // TODO - load test input
      problem_manager.InitCase(
        *eval_hardware,
        org,
        test_id,
        training
      );
      emp_assert(eval_hardware->ValidateThreadState());
    }
  );

  do_program_test_sig.AddAction(
    [this](org_t& org, size_t test_id) {
      emp_assert(eval_hardware->ValidateThreadState());
      // Step the hardware forward to process the input signal
      for (size_t step = 0; step < config.EVAL_CPU_CYCLES_PER_TEST(); ++step) {
        eval_hardware->SingleProcess();
        // Stop early if no active or pending threads
        const size_t num_active_threads = eval_hardware->GetNumActiveThreads();
        const size_t num_pending_threads = eval_hardware->GetNumPendingThreads();
        const bool stop_early = eval_hardware->GetCustomComponent().GetStopEval();
        if (!(num_active_threads || num_pending_threads) || stop_early) {
          break;
        }
      }
    }
  );

  end_program_test_sig.AddAction(
    [this](org_t& org, size_t test_id, bool training) {
      const size_t org_id = org.GetPopID();
      TestResult result = problem_manager.EvaluateOutput(
        *eval_hardware,
        org,
        test_id,
        training
      );
      // Record result on organism phenotype
      auto& phen = org.GetPhenotype();
      // TODO - consolodate phenotype tracking and local world tracking vectors
      phen.test_scores[test_id] = result.score;
      phen.test_passes[test_id] = result.is_correct;
      phen.test_evaluated[test_id] = true;
      phen.aggregate_score += result.score;
      // world performance tracking
      // TODO - assign partial credit for producing an output?
      // org_aggregate_scores[org_id] += result.score;
      org_training_scores[org_id][test_id] = result.score;
      org_training_evaluations[org_id][test_id] = true;
      org_training_coverage[org_id] += (size_t)result.is_correct;
      org_num_training_cases[org_id] += 1;
      // TODO - make sure pop training coverage gets 0'd out at beginning of evaluation step!
      pop_training_coverage[test_id] = pop_training_coverage[test_id] || result.is_correct;
    }
  );

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
    for (size_t test_id = 0; test_id < total_training_cases; ++test_id) {
      fit_fun_set[org_id].emplace_back(
        [this, org_id, test_id]() {
          emp_assert(org_id < org_training_evaluations.size());
          emp_assert(org_id < org_training_scores.size());
          emp_assert(test_id < org_training_evaluations[org_id].size());
          emp_assert(test_id < org_training_scores[org_id].size());
          if (org_training_evaluations[org_id][test_id]) {
            return org_training_scores[org_id][test_id];
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

}

void ProgSynthWorld::SetupEvaluation_Full() {
  std::cout << "Configuring evaluation mode: full" << std::endl;
  emp_assert(total_training_cases > 0);

  // Initialize the test groupings with one group that holds all tests.
  test_groupings.SetSingleGroupMode();
  // Initialize the organism groupings with one group that holds all organisms.
  org_groupings.SetSingleGroupMode();
  emp_assert(test_groupings.GetNumGroups() == org_groupings.GetNumGroups());

  // Configure organism evaluation
  do_org_evaluation_sig.AddAction(
    [this](size_t org_id) {
      emp_assert(org_id < GetSize());
      emp_assert(test_groupings.GetNumGroups() == 1);
      auto& org = GetOrg(org_id);
      begin_program_eval_sig.Trigger(org);
      const auto& test_group = test_groupings.GetGroup(0);
      const auto& test_ids = test_group.GetMembers();
      for (size_t i = 0; i < test_ids.size(); ++i) {
        const size_t test_id = test_ids[i]; // Get test id from group.
        // Handles test input:
        begin_program_test_sig.Trigger(org, test_id, true);
        // Runs the program:
        do_program_test_sig.Trigger(org, test_id);
        // Handles test output evaluation, updates phenotype:
        end_program_test_sig.Trigger(org, test_id, true);
        ++total_test_evaluations;
      }
    }
  );

  // Configure check to see if organism is a solution or not
  check_org_solution = [this](size_t org_id) {
    // TODO - Implement!
    emp_assert(false);
    return false;
  };

}

void ProgSynthWorld::SetupEvaluation_Cohort() {
  // TODO
  emp_assert(false);
}

void ProgSynthWorld::SetupEvaluation_DownSample() {
  // TODO
  emp_assert(false);
}

void ProgSynthWorld::SetupSelection() {
  std::cout << "Configuring parent selection routine" << std::endl;
  // TODO - Can I have different worlds share selection routine setups?
  emp_assert(selector == nullptr);

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

void ProgSynthWorld::SetupSelection_Lexicase() {

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

void ProgSynthWorld::SetupSelection_Tournament() {
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

void ProgSynthWorld::SetupSelection_Truncation() {
  // TODO
  emp_assert(false);
}

void ProgSynthWorld::SetupSelection_None() {
  // TODO
  emp_assert(false);
}

void ProgSynthWorld::SetupSelection_Random() {
  // TODO
  emp_assert(false);
}

void ProgSynthWorld::SetupFitFunEstimator() {
  // Setup the default estimate_test_score functionality
  // TODO - setup configurable fitness estimation
  std::cout << "Configuring fitness function estimator (mode: " << config.EVAL_FIT_EST_MODE() << ")" << std::endl;

  if (config.EVAL_FIT_EST_MODE() == "none") {
    SetupFitFunEstimator_None();
  }
  else if (config.EVAL_FIT_EST_MODE() == "ancestor") {
    SetupFitFunEstimator_Ancestor();
  }
  else if (config.EVAL_FIT_EST_MODE() == "relative") {
    SetupFitFunEstimator_Relative();
  } else {
    std::cout << "Unrecognized EVAL_FIT_EST_MODE: " << config.EVAL_FIT_EST_MODE() << std::endl;
    exit(-1);
  }
}

void ProgSynthWorld::SetupFitFunEstimator_None() {
  // Don't estimate anything
  estimate_test_score = [this](size_t org_id, size_t test_id) {
    return org_training_scores[org_id][test_id];
  };

  // Configure selection routine
  // No estimation, so only use evaluated tests in selection routine
  run_selection_routine = [this]() {
    // Resize parent ids to hold pop_size parents
    selected_parent_ids.resize(config.POP_SIZE(), 0);
    emp_assert(test_groupings.GetNumGroups() == org_groupings.GetNumGroups());
    const size_t num_groups = org_groupings.GetNumGroups();
    // For each grouping, select a number of parents equal to group size
    size_t num_selected = 0;
    for (size_t group_id = 0; group_id < num_groups; ++group_id) {

      auto& org_group = org_groupings.GetGroup(group_id);
      auto& test_group = test_groupings.GetGroup(group_id);
      const size_t n = org_group.GetSize();
      auto& selected = selection_fun(
        n,
        org_group.GetMembers(),
        test_group.GetMembers()
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

}

void ProgSynthWorld::SetupFitFunEstimator_Ancestor() {

  // estimate_test_score
  estimate_test_score = [this](size_t org_id, size_t test_id) {
    emp::Ptr<taxon_t> taxon_ptr = systematics_ptr->GetTaxonAt(org_id);

    auto ancestor = phylo::NearestAncestorWithTraitEval(
      taxon_ptr,
      test_id,
      (size_t)config.EVAL_MAX_PHYLO_SEARCH_DEPTH()
    );

    if (ancestor) {
      auto found_tax = ancestor.value();
      emp_assert(found_tax->GetData().GetTraitsEvaluated()[test_id]);
      return found_tax->GetData().GetPhenotype()[test_id];
    } else {
      return org_training_scores[org_id][test_id];
    }

  };

  // run_selection_routine
  run_selection_routine = [this]() {
    // Resize parent ids to hold pop_size parents
    selected_parent_ids.resize(config.POP_SIZE(), 0);
    emp_assert(test_groupings.GetNumGroups() == org_groupings.GetNumGroups());
    const size_t num_groups = org_groupings.GetNumGroups();
    // For each grouping, select a number of parents equal to group size
    size_t num_selected = 0;
    for (size_t group_id = 0; group_id < num_groups; ++group_id) {
      auto& org_group = org_groupings.GetGroup(group_id);
      const size_t n = org_group.GetSize();
      // Run selection, but use all possible test ids
      auto& selected = selection_fun(
        n,
        org_group.GetMembers(),
        all_training_case_ids
      );
      emp_assert(selected.size() == n);
      emp_assert(n + num_selected <= selected_parent_ids.size());
      std::copy(
        selected.begin(),
        selected.end(),
        selected_parent_ids.begin() + num_selected
      );
      num_selected += n;
    }
  };
}

void ProgSynthWorld::SetupFitFunEstimator_Relative() {

  // estimate_test_score
  estimate_test_score = [this](size_t org_id, size_t test_id) {
    emp::Ptr<taxon_t> taxon_ptr = systematics_ptr->GetTaxonAt(org_id);

    auto ancestor = phylo::NearestRelativeWithTraitEval(
      taxon_ptr,
      test_id,
      (size_t)config.EVAL_MAX_PHYLO_SEARCH_DEPTH()
    );

    if (ancestor) {
      auto found_tax = ancestor.value();
      emp_assert(found_tax->GetData().GetTraitsEvaluated()[test_id]);
      return found_tax->GetData().GetPhenotype()[test_id];
    } else {
      return org_training_scores[org_id][test_id];
    }

  };

  // run_selection_routine
  run_selection_routine = [this]() {
    // Resize parent ids to hold pop_size parents
    selected_parent_ids.resize(config.POP_SIZE(), 0);
    emp_assert(test_groupings.GetNumGroups() == org_groupings.GetNumGroups());
    const size_t num_groups = org_groupings.GetNumGroups();
    // For each grouping, select a number of parents equal to group size
    size_t num_selected = 0;
    for (size_t group_id = 0; group_id < num_groups; ++group_id) {
      auto& org_group = org_groupings.GetGroup(group_id);
      const size_t n = org_group.GetSize();
      // Run selection, but use all possible test ids
      auto& selected = selection_fun(
        n,
        org_group.GetMembers(),
        all_training_case_ids
      );
      emp_assert(selected.size() == n);
      emp_assert(n + num_selected <= selected_parent_ids.size());
      std::copy(
        selected.begin(),
        selected.end(),
        selected_parent_ids.begin() + num_selected
      );
      num_selected += n;
    }

  };

}

void ProgSynthWorld::SetupStoppingCondition() {
  std::cout << "Configuring stopping condition" << std::endl;
  if (config.STOP_MODE() == "generations") {
    SetupStoppingCondition_Generations();
  } else if (config.STOP_MODE() == "evaluations") {
    SetupStoppingCondition_Evaluations();
  } else {
    std::cout << "Unknown STOP_MODE: " << config.STOP_MODE() << std::endl;
    exit(-1);
  }
}

void ProgSynthWorld::SetupStoppingCondition_Generations() {
  stop_run = [this]() {
    return GetUpdate() > config.MAX_GENS();
  };
  is_final_update = [this]() {
    return GetUpdate() >= config.MAX_GENS();
  };
}

void ProgSynthWorld::SetupStoppingCondition_Evaluations() {
  stop_run = [this]() {
    return total_test_evaluations > config.MAX_EVALS();
  };
  is_final_update = [this]() {
    return total_test_evaluations >= config.MAX_EVALS();
  };
}

void ProgSynthWorld::InitializePopulation() {
  std::cout << "Initializing the population." << std::endl;
  // Clear the current population
  Clear();
  // Initialize population according to configuration
  if (config.POP_INIT_MODE() == "random") {
    InitializePopulation_Random();
  } else if (config.POP_INIT_MODE() == "load") {
    InitializePopulation_Load();
  } else {
    std::cout << "Unknown POP_INIT_MODE: " << config.POP_INIT_MODE() << std::endl;
    exit(-1);
  }
}

void ProgSynthWorld::InitializePopulation_Load() {
  // TODO
  emp_assert(false);
}

void ProgSynthWorld::InitializePopulation_Random() {
  for (size_t i = 0; i < config.POP_SIZE(); ++i) {
    Inject(
      {
        sgp::cpu::lfunprg::GenRandLinearFunctionsProgram<hardware_t, TAG_SIZE>(
          *random_ptr,
          inst_lib,
          {config.PRG_MIN_FUNC_CNT(), config.PRG_MAX_FUNC_CNT()},
          FUNC_NUM_TAGS,
          {config.PRG_MIN_FUNC_INST_CNT(), config.PRG_MAX_FUNC_INST_CNT()},
          INST_TAG_CNT,
          INST_ARG_CNT,
          {config.PRG_INST_MIN_ARG_VAL(), config.PRG_INST_MAX_ARG_VAL()}
        )
      }
    );
  }
}

void ProgSynthWorld::SetupPhylogenyTracking() {
  std::cout << "Configure phylogeny tracking" << std::endl;
  emp_assert(systematics_ptr == nullptr);
  // Create new systematics tracker
  systematics_ptr = emp::NewPtr<systematics_t>(
    [](const org_t& org) {
      return org.GetGenome();
    }
  );


  // Add phylogeny snapshot functions
  // Fitness (aggregate score)
  systematics_ptr->AddSnapshotFun(
    [](const taxon_t& taxon) {
      return emp::to_string(taxon.GetData().GetFitness());
    },
    "fitness"
  );

  // Phenotype
  systematics_ptr->AddSnapshotFun(
    [](const taxon_t& taxon) {
      std::stringstream ss;
      utils::PrintVector(ss, taxon.GetData().GetPhenotype(), true);
      return ss.str();
    },
    "phenotype"
  );

  // Genome - TODO
  systematics_ptr->AddSnapshotFun(
    [](const taxon_t& taxon) {
      return "TODO";
    },
    "genome"
  );

  systematics_ptr->AddEvolutionaryDistinctivenessDataNode();
  systematics_ptr->AddPairwiseDistanceDataNode();
  systematics_ptr->AddPhylogeneticDiversityDataNode();

  AddSystematics(systematics_ptr, "genotype");
  SetupSystematicsFile(
    "genotype",
    output_dir + "systematics.csv"
  ).SetTimingRepeat(config.OUTPUT_SUMMARY_DATA_INTERVAL());


}

void ProgSynthWorld::SetupDataCollection() {
  std::cout << "Configure data tracking" << std::endl;
  SetupDataCollection_Phylodiversity();
  SetupDataCollection_Summary();
  SetupDataCollection_Elite();
}

void ProgSynthWorld::SetupDataCollection_Phylodiversity() {
  // Create phylodiversity file
  phylodiversity_file_ptr = emp::NewPtr<emp::DataFile>(
    output_dir + "phylodiversity.csv"
  );

  phylodiversity_file_ptr->AddVar(update, "update", "Generation");
  phylodiversity_file_ptr->AddVar(total_test_evaluations, "evaluations", "Test evaluations so far");
  phylodiversity_file_ptr->AddStats(*systematics_ptr->GetDataNode("evolutionary_distinctiveness") , "genotype_evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
  phylodiversity_file_ptr->AddStats(*systematics_ptr->GetDataNode("pairwise_distance"), "genotype_pairwise_distance", "pairwise distance for a single update", true, true);
  phylodiversity_file_ptr->AddCurrent(*systematics_ptr->GetDataNode("phylogenetic_diversity"), "genotype_current_phylogenetic_diversity", "current phylogenetic_diversity", true, true);
  phylodiversity_file_ptr->PrintHeaderKeys();
}

void ProgSynthWorld::SetupDataCollection_Summary() {
  // Create summary file
  summary_file_ptr = emp::NewPtr<emp::DataFile>(
    output_dir + "summary.csv"
  );

  // TODO
}

void ProgSynthWorld::SetupDataCollection_Elite() {
  // Create elite file
  elite_file_ptr = emp::NewPtr<emp::DataFile>(
    output_dir + "elite.csv"
  );

  // TODO
}

void ProgSynthWorld::SnapshotConfig() {
  emp::DataFile snapshot_file(output_dir + "run_config.csv");
  std::function<std::string(void)> get_param;
  std::function<std::string(void)> get_value;

  snapshot_file.AddFun<std::string>(
    [&get_param]() { return get_param(); },
    "parameter"
  );
  snapshot_file.AddFun<std::string>(
    [&get_value]() { return get_value(); },
    "value"
  );
  snapshot_file.PrintHeaderKeys();

  // Snapshot everything from config file
  for (const auto& entry : config) {
    get_param = [&entry]() { return entry.first; };
    get_value = [&entry]() { return emp::to_string(entry.second->GetValue()); };
    snapshot_file.Update();
  }
  // Snapshot misc. other details
  emp::vector<std::pair<std::string, std::string>> misc_params = {
    std::make_pair("TAG_SIZE", emp::to_string(TAG_SIZE)),
    std::make_pair("FUNC_NUM_TAGS", emp::to_string(FUNC_NUM_TAGS)),
    std::make_pair("INST_TAG_CNT", emp::to_string(INST_TAG_CNT)),
    std::make_pair("INST_ARG_CNT", emp::to_string(INST_ARG_CNT)),
    std::make_pair("sgp_program_type", "LinearFunctionsProgram"),
    std::make_pair("matchbin_metric", "HammingMetric"),
    std::make_pair("matchbin_selector", "RankedSelector"),
    std::make_pair("total_training_cases", emp::to_string(total_training_cases)),
    std::make_pair("total_testing_cases", emp::to_string(problem_manager.GetTestingSetSize()))
  };

  for (const auto& entry : misc_params) {
    get_param = [&entry]() { return entry.first; };
    get_value = [&entry]() { return entry.second; };
    snapshot_file.Update();
  }

}

}