#define CATCH_CONFIG_MAIN

#include "Catch2/single_include/catch2/catch.hpp"

#include "emp/matching/MatchBin.hpp"
#include "emp/bits/BitSet.hpp"

#include "sgp/cpu/mem/BasicMemoryModel.hpp"
#include "sgp/cpu/LinearFunctionsProgramCPU.hpp"

#include "program-synthesis/MutatorLinearFunctionsProgram.hpp"

TEST_CASE("MutatorLinearFunctionsProgram") {
  constexpr size_t TAG_WIDTH = 16;
  constexpr int RANDOM_SEED = 1;
  using mem_model_t = sgp::cpu::mem::BasicMemoryModel;
  using tag_t = emp::BitSet<TAG_WIDTH>;
  using arg_t = int;
  using matchbin_t = emp::MatchBin<
    size_t,
    emp::HammingMetric<TAG_WIDTH>,
    emp::RankedSelector<>,
    emp::NopRegulator
  >;
  using hardware_t = sgp::cpu::LinearFunctionsProgramCPU<
    mem_model_t,
    arg_t,
    matchbin_t
  >;
  using inst_t = typename hardware_t::inst_t;
  using inst_lib_t = typename hardware_t::inst_lib_t;
  using program_t = typename hardware_t::program_t;
  using mutator_t = psynth::MutatorLinearFunctionsProgram<
    hardware_t,
    tag_t,
    arg_t
  >;

  inst_lib_t inst_lib;

  inst_lib.AddInst("Nop-A", [](hardware_t & hw, const inst_t & inst) { ; }, "No operation!");
  inst_lib.AddInst("Nop-B", [](hardware_t & hw, const inst_t & inst) { ; }, "No operation!");
  inst_lib.AddInst("Nop-C", [](hardware_t & hw, const inst_t & inst) { ; }, "No operation!");
  emp::Random random(RANDOM_SEED);
  mutator_t mutator(inst_lib);
  emp::Range<int> ARG_VAL_RANGE = {0, 15};
  emp::Range<size_t> FUNC_LEN_RANGE = {0, 128};
  emp::Range<size_t> FUNC_CNT_RANGE = {1, 32};
  size_t MAX_TOTAL_LEN = 1024;
  size_t NUM_FUNC_TAGS = 1;
  size_t NUM_INST_TAGS = 1;
  size_t NUM_INST_ARGS = 3;
  mutator.SetProgFunctionCntRange(FUNC_CNT_RANGE);
  mutator.SetProgFunctionInstCntRange(FUNC_LEN_RANGE);
  mutator.SetProgInstArgValueRange(ARG_VAL_RANGE);
  mutator.SetTotalInstLimit(MAX_TOTAL_LEN);
  mutator.SetFuncNumTags(NUM_FUNC_TAGS);
  mutator.SetInstNumTags(NUM_INST_TAGS);
  mutator.SetInstNumArgs(NUM_INST_ARGS);
  // Test 0 mutation rate on all functions.
  mutator.SetRateInstArgSub(0.0);
  mutator.SetRateInstTagBF(0.0);
  mutator.SetRateInstSub(0.0);
  mutator.SetRateInstIns(0.0);
  mutator.SetRateInstDel(0.0);
  mutator.SetRateSeqSlip(0.0);
  mutator.SetRateFuncDup(0.0);
  mutator.SetRateFuncDel(0.0);
  mutator.SetRateFuncTagBF(0.0);
  program_t nop_prog;
  size_t num_muts = 0;
  for (size_t f = 0; f < 3; ++f) {
    nop_prog.PushFunction(tag_t());
    for (size_t i = 0; i < 8; ++i) {
      nop_prog.PushInst(inst_lib, "Nop-A", {0, 0, 0}, {tag_t()});
    }
  }
  program_t copy_prog(nop_prog);
  num_muts = mutator.ApplyInstSubs(random, nop_prog);
  REQUIRE(num_muts == 0);
  REQUIRE(copy_prog == nop_prog);
  REQUIRE(mutator.VerifyProgram(nop_prog));
  num_muts = mutator.ApplyInstInDels(random, nop_prog);
  REQUIRE(num_muts == 0);
  REQUIRE(copy_prog == nop_prog);
  REQUIRE(mutator.VerifyProgram(nop_prog));
  num_muts = mutator.ApplySeqSlips(random, nop_prog);
  REQUIRE(num_muts == 0);
  REQUIRE(copy_prog == nop_prog);
  REQUIRE(mutator.VerifyProgram(nop_prog));
  num_muts = mutator.ApplyFuncDup(random, nop_prog);
  REQUIRE(num_muts == 0);
  REQUIRE(copy_prog == nop_prog);
  REQUIRE(mutator.VerifyProgram(nop_prog));
  num_muts = mutator.ApplyFuncDel(random, nop_prog);
  REQUIRE(num_muts == 0);
  REQUIRE(copy_prog == nop_prog);
  REQUIRE(mutator.VerifyProgram(nop_prog));
  num_muts = mutator.ApplyFuncTagBF(random, nop_prog);
  REQUIRE(num_muts == 0);
  REQUIRE(copy_prog == nop_prog);
  REQUIRE(mutator.VerifyProgram(nop_prog));
  // Check function duplications.
  mutator.SetRateFuncDup(1.0);
  size_t orig_f_cnt = nop_prog.GetSize();
  mutator.ApplyAll(random, nop_prog);
  REQUIRE(nop_prog.GetSize() == 2*orig_f_cnt);
  REQUIRE(mutator.VerifyProgram(nop_prog));
  // Check function deletions.
  mutator.SetRateFuncDel(1.0);
  mutator.SetRateFuncDup(0.0);
  mutator.ApplyAll(random, nop_prog);
  REQUIRE(nop_prog.GetSize() == FUNC_CNT_RANGE.GetLower());
  REQUIRE(mutator.VerifyProgram(nop_prog));
  // Generate many random programs, apply mutations, check constraints.
  mutator.SetRateInstArgSub(0.25);
  mutator.SetRateInstTagBF(0.25);
  mutator.SetRateInstSub(0.25);
  mutator.SetRateInstIns(0.25);
  mutator.SetRateInstDel(0.25);
  mutator.SetRateSeqSlip(0.25);
  mutator.SetRateFuncDup(0.25);
  mutator.SetRateFuncDel(0.25);
  mutator.SetRateFuncTagBF(0.25);
  for (size_t i = 0; i < 1000; ++i) {
    program_t prog(
      sgp::cpu::lfunprg::GenRandLinearFunctionsProgram<hardware_t, TAG_WIDTH>(
        random,
        inst_lib,
        {1,16},
        NUM_FUNC_TAGS,
        {0,64},
        NUM_INST_TAGS,
        NUM_INST_ARGS,
        ARG_VAL_RANGE
      )
    );
    REQUIRE(mutator.VerifyProgram(prog));
    for (size_t m = 0; m < 100; ++m) {
      mutator.ApplyAll(random, prog);
      REQUIRE(mutator.VerifyProgram(prog));
    }
  }
}