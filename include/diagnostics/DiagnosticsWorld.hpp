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

protected:
  const config_t& config;

  emp::Ptr<BaseDiagnostic> base_diagnostic=nullptr; ///< Base-layer diagnostic to use to translate genomes to phenotypes
  MultiValleyCrossingDiagnostic valley_diagnostic;  ///< If in use, it's layered on top of another diagnostic (base_diagnostic).

  std::function<void(org_t&)> evaluate_org_fun;
  std::function<void(const genome_t&, phenotype_t&)> translate_genome_fun;

  void Setup();
  void SetupDiagnostic();
  void SetupSelection();
  void SetupMutator();
  void SetupDataCollection();

  template<typename DIAG_PROB>
  void SetupDiagnosticHelper();

  void SetupDiagnostic_ExploitationRate();
  void SetupDiagnostic_OrderedExploitation();
  void SetupDiagnostic_ContradictoryObjectives();
  void SetupDiagnostic_MultipathExploration();

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

  // Configure organism evaluation
  evaluate_org_fun = [this](org_t& org) {
    // Translate organism genome
    org.TranslateGenome(translate_genome_fun);
    // Calculate optimal traits
    org.CalcOptimalTraits(config.TARGET(), config.ACCURACY());
  };
}



}