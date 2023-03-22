#define CATCH_CONFIG_MAIN

#include "Catch2/single_include/catch2/catch.hpp"

#include "emp/matching/MatchBin.hpp"
#include "emp/bits/BitSet.hpp"

#include "sgp/cpu/mem/BasicMemoryModel.hpp"
#include "sgp/cpu/LinearFunctionsProgramCPU.hpp"
#include "sgp/cpu/lfunprg/LinearFunctionsProgram.hpp"
#include "sgp/inst/lfpbm/InstructionAdder.hpp"

#include "program-synthesis/program_utils.hpp"

TEST_CASE("PrintProgram") {
  constexpr size_t TAG_WIDTH = 16;
  constexpr int RANDOM_SEED = 2;
  using mem_model_t = sgp::cpu::mem::BasicMemoryModel;
  using arg_t = int;
  using matchbin_t = emp::MatchBin<
    size_t,
    emp::HammingMetric<TAG_WIDTH>,
    emp::RankedSelector<>,
    emp::AdditiveCountdownRegulator<>
  >;
  using hardware_t = sgp::cpu::LinearFunctionsProgramCPU<
    mem_model_t,
    arg_t,
    matchbin_t
  >;
  using inst_lib_t = typename hardware_t::inst_lib_t;
  using program_t = typename hardware_t::program_t;
  // using function_t = typename program_t::function_t;
  using tag_t = emp::BitSet<TAG_WIDTH>;

  emp::Random random(RANDOM_SEED);

  inst_lib_t inst_lib;
  inst_lib.Clear();
  sgp::inst::lfpbm::InstructionAdder<hardware_t> inst_adder;
  // Add instructions (except Fork and Terminate)
  inst_adder.AddAllDefaultInstructions(
    inst_lib,
    {"Fork", "Terminate"}
  );

  // Build a program
  program_t program;
  program.PushFunction(tag_t(random));
  program.PushInst(
    inst_lib,
    "Nop",
    {0, 1, 2},
    {tag_t(random)}
  );
  program.PushInst(
    inst_lib,
    "Nop",
    {3, 4, 5},
    {tag_t(random)}
  );
  program.PushInst(
    inst_lib,
    "Nop",
    {5, 6, 7},
    {tag_t(random)}
  );
  program.PushFunction(emp::vector<tag_t>{tag_t(random), tag_t(random), tag_t(random)});
  program.PushInst(
    inst_lib,
    "Add",
    {0, 1, 2},
    {tag_t(random)}
  );
  program.PushFunction(tag_t(random));
  program.PushInst(
    inst_lib,
    "Mult",
    {0, 1, 2},
    {tag_t(random)}
  );

  SECTION("PrintProgramJSON") {
    psynth::PrintProgramJSON(
      std::cout,
      program,
      inst_lib
    );
    std::cout << std::endl;
  }

  SECTION("Print") {
    program.Print(
      std::cout,
      inst_lib
    );
  }
}

TEST_CASE("LoadProgram") {
  constexpr size_t TAG_WIDTH = 16;
  // constexpr int RANDOM_SEED = 2;
  using mem_model_t = sgp::cpu::mem::BasicMemoryModel;
  using arg_t = int;
  using matchbin_t = emp::MatchBin<
    size_t,
    emp::HammingMetric<TAG_WIDTH>,
    emp::RankedSelector<>,
    emp::AdditiveCountdownRegulator<>
  >;
  using hardware_t = sgp::cpu::LinearFunctionsProgramCPU<
    mem_model_t,
    arg_t,
    matchbin_t
  >;
  using inst_lib_t = typename hardware_t::inst_lib_t;
  using program_t = typename hardware_t::program_t;
  // using function_t = typename program_t::function_t;
  // using tag_t = emp::BitSet<TAG_WIDTH>;

  // emp::Random random(RANDOM_SEED);

  inst_lib_t inst_lib;
  inst_lib.Clear();
  sgp::inst::lfpbm::InstructionAdder<hardware_t> inst_adder;
  // Add instructions (except Fork and Terminate)
  inst_adder.AddAllDefaultInstructions(
    inst_lib,
    {"Fork", "Terminate"}
  );

  std::ifstream prg_fstream("sgp_program.sgp");
  REQUIRE(prg_fstream.is_open());

  program_t program = psynth::LoadLinearFunctionsProgram_PrintFormat<inst_lib_t, TAG_WIDTH>(
    prg_fstream,
    inst_lib
  );

  program.Print(
    std::cout,
    inst_lib
  );
}