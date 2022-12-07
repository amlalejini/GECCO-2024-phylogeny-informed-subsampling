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
  using config_t = DiagnosticsConfig;
protected:
  const config_t& config;

  void Setup();
  void SetupSelection();
  void SetupMutator();
  void SetupDataCollection();

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
  }

  void RunStep();
  void Run();

};

void DiagnosticsWorld::Setup() {
  std::cout << "--- Setting up DiagnosticsWorld ---" << std::endl;
}



}