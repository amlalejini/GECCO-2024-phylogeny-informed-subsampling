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

#include "sgp/cpu/linprg/LinearProgram.hpp"
#include "sgp/cpu/LinearProgramCPU.hpp"
#include "sgp/cpu/mem/BasicMemoryModel.hpp"
#include "sgp/inst/InstructionLibrary.hpp"
#include "sgp/EventLibrary.hpp"
#include "sgp/inst/lpbm/InstructionAdder.hpp"

#include "../phylogeny/phylogeny_utils.hpp"

#include "ProgSynthConfig.hpp"
#include "ProgSynthOrg.hpp"
#include "Event.hpp"
#include "ProblemManager.hpp"
#include "ProgSynthHardware.hpp"

namespace psynth {

// TODO - move these into internal namespace
constexpr size_t TAG_SIZE = 32;
using tag_t = emp::BitSet<TAG_SIZE>;
using inst_arg_t = int;
using program_t = sgp::cpu::linprg::LinearProgram<tag_t, inst_arg_t>;
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
  using hardware_t = sgp::cpu::LinearProgramCPU<
    hw_memory_model_t,
    inst_arg_t,
    hw_matchbin_t,
    ProgSynthHardwareComponent
  >;
  using inst_lib_t = sgp::inst::InstructionLibrary<hardware_t, inst_t>;
  using event_lib_t = sgp::EventLibrary<hardware_t>;
  using base_event_t = typename event_lib_t::event_t;
  using systematics_t = emp::Systematics<
    org_t,
    genome_t,
    phylo::phenotype_info
  >;
  using taxon_info_t = phylo::phenotype_info;
  using taxon_t = typename systematics_t::taxon_t;

  using config_t = ProgSynthConfig;
  using selection_fun_t = std::function<
    emp::vector<size_t>&(
      size_t,
      const emp::vector<size_t>&,
      const emp::vector<size_t>&
    )
  >;

  struct Grouping {
    size_t group_id = 0;
    size_t group_size = 0;
    emp::vector<size_t> member_ids;

    void Resize(size_t size, size_t value = 0) {
      group_size = size;
      member_ids.resize(group_size, value);
    }

    size_t GetSize() const { return member_ids.size(); }
  };

  struct SelectedStatistics {
    // TODO
  };

protected:
  const config_t& config;

  size_t total_test_evaluations = 0; ///< Tracks the total number of "test case" evaluations across all organisms since the beginning of the run.

  emp::Ptr<hardware_t> eval_hardware;
  inst_lib_t inst_lib;
  event_lib_t event_lib;

  ProblemManager problem_manager;

  size_t event_id_input_sig = 0;

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

  // Setup the population structure
  SetPopStruct_Mixed(true);

  // Configure world to set organism ID on placement
  OnPlacement(
    [this](size_t pos) {
      auto& org = GetOrg(pos);
      org.SetPopID(pos);
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
  // Reset instruction library
  inst_lib.Clear();
  sgp::inst::lpbm::InstructionAdder<hardware_t> inst_adder;
  inst_adder.AddAllDefaultInstructions(inst_lib);
  // TODO - add problem-specific instructions

  // TODO - snapshot instruction set
}

void ProgSynthWorld::SetupEventLibrary() {
  event_lib.Clear();
  // TODO - configure event library
  // event_id_input_sig = event_lib.AddEvent(
  //   "InputSignal",
  //   [this](hardware_t& hw, const base_event_t& e) {

  //   }
  // );
}

void ProgSynthWorld::SetupVirtualHardware() {
  if (eval_hardware == nullptr) {
    eval_hardware = emp::NewPtr<hardware_t>(*random_ptr, inst_lib, event_lib);
  }
  // Configure the SGP CPU
  // TODO - Implement ability to run hardware in parallel
  eval_hardware->Reset();
  eval_hardware->SetActiveThreadLimit(config.MAX_ACTIVE_THREAD_CNT());
  eval_hardware->SetThreadCapacity(config.MAX_THREAD_CAPACITY());
  emp_assert(eval_hardware->ValidateThreadState());
}

}