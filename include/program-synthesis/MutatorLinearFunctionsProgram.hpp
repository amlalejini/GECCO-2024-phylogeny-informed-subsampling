/* MutatorLinearFunctionsProgram.hpp

Implemenation adapted from: https://github.com/amlalejini/Tag-based-Genetic-Regulation-for-LinearGP/

*/
#pragma once

#include <unordered_map>
#include "emp/bits/BitSet.hpp"
#include "emp/math/Random.hpp"
#include "emp/math/random_utils.hpp"
#include "emp/math/Range.hpp"

#include "sgp/cpu/lfunprg/LinearFunctionsProgram.hpp"

namespace psynth {

template<
  typename HARDWARE_T,
  typename TAG_T,
  typename ARGUMENT_T
>
class MutatorLinearFunctionsProgram {
public:

  MutatorLinearFunctionsProgram() {
    emp_assert(false, "Not implemented!");
  }

};

/// Mutator for SignalGP programs that use bitstrings for tags and ints as instruction arguments.
// the slow add-on of mutation types has caused this class to bloat. Future projects should re-design
// the mutator (to make managing many mutation types simpler).
template<typename HARDWARE_T, size_t TAG_W>
class MutatorLinearFunctionsProgram<HARDWARE_T, emp::BitSet<TAG_W>, int> {
public:
  using tag_t = emp::BitSet<TAG_W>;
  using arg_t = int;
  using program_t = sgp::cpu::lfunprg::LinearFunctionsProgram<tag_t, arg_t>;
  using function_t = typename program_t::function_t;
  using inst_t = typename program_t::inst_t;
  using hardware_t = HARDWARE_T;
  using inst_lib_t = typename HARDWARE_T::inst_lib_t;

  enum class MUTATION_TYPES {
    INST_ARG_SUB = 0,
    INST_TAG_BIT_FLIP,
    INST_TAG_BIT_SEQ_RANDOMIZATION,
    INST_SUB,
    INST_INS,
    INST_DEL,
    SEQ_SLIP_DUP,
    SEQ_SLIP_DEL,
    FUNC_DUP,
    FUNC_DEL,
    FUNC_TAG_BIT_FLIP,
    FUNC_TAG_BIT_SEQ_RANDOMIZATION
  };

protected:
  // Need an instruction library to know valid
  inst_lib_t& inst_lib;

  // Program constraints.
  emp::Range<size_t> prog_func_cnt_range = {0, 4}; ///< min/max # functions in a program
  emp::Range<size_t> prog_func_inst_range = {0, 8}; ///< min/max # of instructions a function can have.
  emp::Range<int> prog_inst_arg_val_range = {0, 15};
  size_t prog_func_num_tags=1;  ///< Number of tags associated with each function.
  size_t prog_total_inst=4*8;   ///< How many total instructions is a program allowed to have?
  size_t prog_inst_num_args=3;  ///< How many arguments should each instruction have?
  size_t prog_inst_num_tags=1;  ///< How many tags should each instruction have?

  // Rates
  double rate_inst_arg_sub=0.0;
  double rate_inst_tag_bit_flips=0.0;
  double rate_inst_sub=0.0;
  double rate_inst_ins=0.0;
  double rate_inst_del=0.0;
  double rate_seq_slip=0.0;
  double rate_func_dup=0.0;
  double rate_func_del=0.0;
  double rate_func_tag_bit_flips=0.0;
  double rate_func_tag_single_bit_flip=0.0; ///< Per-tag
  double rate_inst_tag_single_bit_flip=0.0; ///< Probability of a *single* bit flip (in random location), applied per-tag.
  double rate_func_tag_seq_rand=0.0;  ///< Per-tag
  double rate_inst_tag_seq_rand=0.0;  ///< Per-tag

  std::unordered_map<MUTATION_TYPES, int> last_mutation_tracker;

  // emp::vector<std::function<size_t(emp::Random &, program_t &)>> active_mutations;

public:
  MutatorLinearFunctionsProgram(inst_lib_t& ilib)
    : inst_lib(ilib)
  {
    ResetLastMutationTracker();
  }

  void ResetLastMutationTracker() {
    last_mutation_tracker[MUTATION_TYPES::INST_ARG_SUB] = 0;
    last_mutation_tracker[MUTATION_TYPES::INST_TAG_BIT_FLIP] = 0;
    last_mutation_tracker[MUTATION_TYPES::INST_SUB] = 0;
    last_mutation_tracker[MUTATION_TYPES::INST_INS] = 0;
    last_mutation_tracker[MUTATION_TYPES::INST_DEL] = 0;
    last_mutation_tracker[MUTATION_TYPES::SEQ_SLIP_DUP] = 0;
    last_mutation_tracker[MUTATION_TYPES::SEQ_SLIP_DEL] = 0;
    last_mutation_tracker[MUTATION_TYPES::FUNC_DUP] = 0;
    last_mutation_tracker[MUTATION_TYPES::FUNC_DEL] = 0;
    last_mutation_tracker[MUTATION_TYPES::FUNC_TAG_BIT_FLIP] = 0;
  }

  const std::unordered_map<MUTATION_TYPES, int>& GetLastMutations() const {
    return last_mutation_tracker;
  }
  std::unordered_map<MUTATION_TYPES, int>& GetLastMutations() {
    return last_mutation_tracker;
  }

  void SetProgFunctionCntRange(const emp::Range<size_t>& val) { prog_func_cnt_range = val; }
  void SetProgFunctionInstCntRange(const emp::Range<size_t>& val) { prog_func_inst_range = val; }
  void SetProgInstArgValueRange(const emp::Range<int>& val) { prog_inst_arg_val_range = val; }
  void SetTotalInstLimit(size_t val) { prog_total_inst = val; }
  void SetFuncNumTags(size_t val) { prog_func_num_tags = val; }
  void SetInstNumArgs(size_t val) { prog_inst_num_args = val; }
  void SetInstNumTags(size_t val) { prog_inst_num_tags = val; }
  void SetRateInstArgSub(double val) { rate_inst_arg_sub = val; }
  void SetRateInstTagBF(double val) { rate_inst_tag_bit_flips = val; }
  void SetRateInstSub(double val) { rate_inst_sub = val; }
  void SetRateInstIns(double val) { rate_inst_ins = val; }
  void SetRateInstDel(double val) { rate_inst_del = val; }
  void SetRateSeqSlip(double val) { rate_seq_slip = val; }
  void SetRateFuncDup(double val) { rate_func_dup = val; }
  void SetRateFuncDel(double val) { rate_func_del = val; }
  void SetRateFuncTagBF(double val) { rate_func_tag_bit_flips = val; }

  void SetRateFuncTagSingleBF(double val) { rate_func_tag_single_bit_flip = val; }
  void SetRateInstTagSingleBF(double val) { rate_inst_tag_single_bit_flip = val; }
  void SetRateFuncTagSeqRand(double val) { rate_func_tag_seq_rand = val; }
  void SetRateInstTagSeqRand(double val) { rate_inst_tag_seq_rand = val; }

  const emp::Range<size_t>& GetProgFunctionCntRange() const { return prog_func_cnt_range; }
  const emp::Range<size_t>& GetProgFunctionInstCntRange() const { return prog_func_inst_range; }
  const emp::Range<int>& GetProgInstArgValueRange() const { return prog_inst_arg_val_range; }
  size_t GetTotalInstLimit() const { return prog_total_inst; }
  size_t GetFuncNumTags() const { return prog_func_num_tags; }
  size_t GetInstNumArgs() const { return prog_inst_num_args; }
  size_t GetInstNumTags() const { return prog_inst_num_tags; }
  double GetRateInstArgSub() const { return rate_inst_arg_sub; }
  double GetRateInstTagBF() const { return rate_inst_tag_bit_flips; }
  double GetRateInstSub() const { return rate_inst_sub; }
  double GetRateInstIns() const { return rate_inst_ins; }
  double GetRateInstDel() const { return rate_inst_del; }
  double GetRateSeqSlip() const { return rate_seq_slip; }
  double GetRateFuncDup() const { return rate_func_dup; }
  double GetRateFuncDel() const { return rate_func_del; }
  double GetRateFuncTagBF() const { return rate_func_tag_bit_flips; }

  double GetRateFuncTagSingleBF() const { return rate_func_tag_single_bit_flip; }
  double GetRateInstTagSingleBF() const { return rate_inst_tag_single_bit_flip; }
  double GetRateFuncTagSeqRand() const { return rate_func_tag_seq_rand; }
  double GetRateInstTagSeqRand() const { return rate_inst_tag_seq_rand; }

  /// Applies all mutation operators at current rates.
  /// Will reset mutation tracking before applying mutations.
  size_t ApplyAll(emp::Random& rnd, program_t& program) {
    size_t mut_cnt = 0;
    mut_cnt += ApplyInstSubs(rnd, program);
    mut_cnt += ApplyInstInDels(rnd, program);
    mut_cnt += ApplySeqSlips(rnd, program);
    mut_cnt += ApplyFuncDup(rnd, program);
    mut_cnt += ApplyFuncDel(rnd, program);
    mut_cnt += ApplyFuncTagBF(rnd, program);
    return mut_cnt;
  }

  /// Apply bit flips to tag @ per-bit rate.
  size_t ApplyTagBitFlipsPerBit(emp::Random& rnd, tag_t& tag, double rate) {
    size_t mut_cnt = 0;
    for (size_t k = 0; k < tag.GetSize(); ++k) {
      if (rnd.P(rate)) {
        tag.Toggle(k);
        ++mut_cnt;
      }
    }
    return mut_cnt;
  }

  /// Apply specified number of tag bit flips.
  size_t ApplyTagBitFlipsFixed(emp::Random& rnd, tag_t& tag, size_t num_flips) {
    tag.FlipRandomCount(rnd, num_flips);
    return num_flips;
  }

  /// Pick two random locations in the tag, randomize everything in between.
  size_t ApplyTagSeqRandomization(emp::Random& rnd, tag_t& tag) {
    const size_t a = rnd.GetUInt(tag.GetSize());
    const size_t b = rnd.GetUInt(tag.GetSize());
    const size_t pos_1 = emp::Min(a, b);
    const size_t pos_2 = emp::Max(a, b);
    emp_assert(pos_1 < tag.size());
    emp_assert(pos_2 < tag.size());
    emp_assert(pos_1 <= pos_2);
    // Mutate positions [pos_1:pos_2] (inclusive)
    // NOTE - we could speed this up with bit magic
    for (size_t pos = pos_1; pos <= pos_2; ++pos) {
      tag.Set(pos, rnd.P(0.5));
    }
    return (pos_2 - pos_1) + 1;
  }

  /// Apply instruction substitutions (operator, argument, tag).
  size_t ApplyInstSubs(emp::Random& rnd, program_t& program) {
    size_t mut_cnt = 0;
    for (size_t fID = 0; fID < program.GetSize(); ++fID) {
      for (size_t iID = 0; iID < program[fID].GetSize(); ++iID) {
        inst_t & inst = program[fID][iID];

        // Mutate instruction tag(s).
        size_t tag_bf_cnt = 0;
        for (tag_t & tag : inst.GetTags()) {
          // Apply per-bit substitution mutations
          if (rate_inst_tag_bit_flips > 0.0) {
            tag_bf_cnt += ApplyTagBitFlipsPerBit(rnd, tag, rate_inst_tag_bit_flips);
          }
          // Apply per-tag substitution mutations
          if (rate_inst_tag_single_bit_flip > 0.0 && rnd.P(rate_inst_tag_single_bit_flip)) {
            tag_bf_cnt += ApplyTagBitFlipsFixed(rnd, tag, 1);
          }
          // Apply per-tag sequence randomization mutations
          if (rate_inst_tag_seq_rand > 0.0 && rnd.P(rate_inst_tag_seq_rand)) {
            ApplyTagSeqRandomization(rnd, tag);
            ++mut_cnt;  // Count this as only one mutation
            ++last_mutation_tracker[MUTATION_TYPES::INST_TAG_BIT_SEQ_RANDOMIZATION];
          }
        }
        mut_cnt += tag_bf_cnt;
        last_mutation_tracker[MUTATION_TYPES::INST_TAG_BIT_FLIP] += tag_bf_cnt;

        // Mutate instruction operation.
        if (rnd.P(rate_inst_sub)) {
          inst.id = rnd.GetUInt(inst_lib.GetSize());
          ++last_mutation_tracker[MUTATION_TYPES::INST_SUB];
          ++mut_cnt;
        }

        // Mutate instruction arguments.
        for (size_t k = 0; k < inst.GetArgs().size(); ++k) {
          if (rnd.P(rate_inst_arg_sub)) {
            inst.GetArgs()[k] = rnd.GetInt(
              prog_inst_arg_val_range.GetLower(),
              prog_inst_arg_val_range.GetUpper()+1
            );
            ++mut_cnt;
            ++last_mutation_tracker[MUTATION_TYPES::INST_ARG_SUB];
          }
        }
      }
    }
    return mut_cnt;
  }

  /// Apply single-instruction insertions and deletions.
  size_t ApplyInstInDels(emp::Random & rnd, program_t & program) {
    size_t mut_cnt = 0;
    size_t expected_prog_len = program.GetInstCount();
    // Perform single-instruction insertion/deletions.
    for (size_t fID = 0; fID < program.GetSize(); ++fID) {
      function_t new_function(program[fID].GetTags()); // Copy over tags
      size_t expected_func_len = program[fID].GetSize();
      // Compute number and location of insertions.
      const uint32_t num_ins = rnd.GetRandBinomial(program[fID].GetSize(), rate_inst_ins);
      emp::vector<size_t> ins_locs;
      if (num_ins > 0) {
        if (program[fID].GetSize()) {
          ins_locs = emp::RandomUIntVector(rnd, num_ins, 0, program[fID].GetSize());
          // std::sort(ins_locs.rbegin(), ins_locs.rend());
          std::sort(ins_locs.begin(), ins_locs.end(), std::greater<size_t>());
        } else {
          ins_locs.resize(num_ins, 0); // If function len is 0, just put all insertions at beginning.
        }
      }
      size_t read_head = 0;
      while (read_head < program[fID].GetSize()) {
        // Should we insert?
        if (ins_locs.size() > 0) {
          if (read_head >= ins_locs.back() &&
              expected_func_len < prog_func_inst_range.GetUpper() &&
              expected_prog_len < prog_total_inst)
          {
            // Insert a new random instruction.
            new_function.PushInst(
              sgp::cpu::linprg::GenRandInst<hardware_t, TAG_W>(
                rnd,
                inst_lib,
                prog_inst_num_tags,
                prog_inst_num_args,
                prog_inst_arg_val_range
              )
            );
            ++mut_cnt;
            ++last_mutation_tracker[MUTATION_TYPES::INST_INS];
            ++expected_func_len;
            ++expected_prog_len;
            ins_locs.pop_back();
            continue;
          }
        }
        // Should we delete this instruction?
        if (rnd.P(rate_inst_del) && expected_func_len > prog_func_inst_range.GetLower()) {
          ++mut_cnt;
          ++last_mutation_tracker[MUTATION_TYPES::INST_DEL];
          --expected_func_len;
          --expected_prog_len;
        } else {
          new_function.PushInst(program[fID][read_head]);
        }
        ++read_head;
      }
      program[fID] = new_function;
    }
    return mut_cnt;
  }

  /// Apply slip mutations to program sequences (per-function instruction sequence).
  size_t ApplySeqSlips(emp::Random& rnd, program_t& program) {
    size_t mut_cnt =0;
    size_t expected_prog_len = program.GetInstCount();
    // Perform per-function slip mutations.
    for (size_t fID = 0; fID < program.GetSize(); ++fID) {
      if (!rnd.P(rate_seq_slip) || program[fID].GetSize() == 0) continue; // don't do it here
      size_t begin = rnd.GetUInt(program[fID].GetSize());
      size_t end = rnd.GetUInt(program[fID].GetSize());
      const bool dup = begin < end;
      const bool del = begin > end;
      const int dup_size = (int)end - (int)begin;
      const int del_size = (int)begin - (int)end;
      if (
        dup &&
        (expected_prog_len + (size_t)dup_size <= prog_total_inst) &&
        (program[fID].GetSize() + (size_t)dup_size <= prog_func_inst_range.GetUpper())
      ) {
        // Duplicate begin:end
        const size_t new_size = program[fID].GetSize() + (size_t)dup_size;
        function_t new_function(program[fID].GetTags());
        for (size_t i = 0; i < new_size; ++i) {
          if (i < end) new_function.PushInst(program[fID][i]);
          else new_function.PushInst(program[fID][i - (size_t)dup_size]);
        }
        program[fID] = new_function;
        ++mut_cnt;
        ++last_mutation_tracker[MUTATION_TYPES::SEQ_SLIP_DUP];
        expected_prog_len += (size_t)dup_size;
      } else if (del && (program[fID].GetSize() - (size_t)del_size) >= prog_func_inst_range.GetLower()) {
        // Delete end:begin
        function_t new_function(program[fID].GetTags());
        for (size_t i = 0; i < end; ++i)
          new_function.PushInst(program[fID][i]);
        for (size_t i = begin; i < program[fID].GetSize(); ++i)
          new_function.PushInst(program[fID][i]);
        program[fID] = new_function;
        ++mut_cnt;
        ++last_mutation_tracker[MUTATION_TYPES::SEQ_SLIP_DEL];
        expected_prog_len -= (size_t)del_size;
      }
    }
   return mut_cnt;
  }

  /// Apply function duplications to program (per-function).
  size_t ApplyFuncDup(emp::Random& rnd, program_t& program) {
    size_t mut_cnt = 0;
    size_t expected_prog_len = program.GetInstCount();
    // Perform function duplications!
    size_t orig_func_wall = program.GetSize();
    for (size_t fID = 0; fID < orig_func_wall; ++fID) {
      // Should we duplicate this function?
      if (fID < orig_func_wall &&
          rnd.P(rate_func_dup) &&
          (program.GetSize() < prog_func_cnt_range.GetUpper()) &&
          (expected_prog_len + program[fID].GetSize() <= prog_total_inst))
      {
        // Duplicate!
        program.PushFunction(program[fID]);
        expected_prog_len += program[fID].GetSize();
        ++mut_cnt;
        ++last_mutation_tracker[MUTATION_TYPES::FUNC_DUP];
      }
    }
    return mut_cnt;
  }

  /// Apply function deletions to program (per-function).
  size_t ApplyFuncDel(emp::Random & rnd, program_t & program) {
    size_t mut_cnt = 0;
    size_t expected_prog_len = program.GetInstCount();
    // Perform function deletions!
    for (int fID = 0; fID < (int)program.GetSize(); ++fID) {
      // Should we delete this function?
      if (rnd.P(rate_func_del) &&
          program.GetSize() > prog_func_cnt_range.GetLower())
      {
        expected_prog_len -= program[(size_t)fID].GetSize();
        program[(size_t)fID] = program[program.GetSize() - 1];
        program.PopFunction();
        ++mut_cnt;
        ++last_mutation_tracker[MUTATION_TYPES::FUNC_DEL];
        fID -= 1;
        continue;
      }
    }
    return mut_cnt;
  }

  /// Apply function tag bit-flip mutations.
  size_t ApplyFuncTagBF(emp::Random& rnd, program_t& program) {
    size_t mut_cnt = 0;
    // Perform function tag mutations!
    for (size_t fID = 0; fID < program.GetSize(); ++fID) {
      size_t tag_bfs = 0;
      for (tag_t & tag : program[fID].GetTags()) {
        // Apply per-bit substitution mutations
        if (rate_func_tag_bit_flips > 0.0) {
          tag_bfs += ApplyTagBitFlipsPerBit(rnd, tag, rate_func_tag_bit_flips);
        }
        // Apply per-tag single-bit substitutions
        if (rate_func_tag_single_bit_flip > 0.0 && rnd.P(rate_func_tag_single_bit_flip)) {
          tag_bfs += ApplyTagBitFlipsFixed(rnd, tag, 1);
        }
        // Apply per-tag sequence randomization substitutions
        if (rate_func_tag_seq_rand > 0.0 && rnd.P(rate_func_tag_seq_rand)) {
          ApplyTagSeqRandomization(rnd, tag);
          ++mut_cnt;
          ++last_mutation_tracker[MUTATION_TYPES::FUNC_TAG_BIT_SEQ_RANDOMIZATION];
        }
      }
      mut_cnt += tag_bfs;
      last_mutation_tracker[MUTATION_TYPES::FUNC_TAG_BIT_FLIP] += tag_bfs;
    }
    return mut_cnt;
  }

  /// Verify that the given program (prog) is within the constraints associated with this SignalGPMutator object.
  /// Useful for mutator testing.
  bool VerifyProgram(program_t & prog) {
    if (prog.GetInstCount() > prog_total_inst) { return false; }
    if (!prog_func_cnt_range.Valid(prog.GetSize())) { return false; }
    for (size_t fID = 0; fID < prog.GetSize(); ++fID) {
      if (!prog_func_inst_range.Valid(prog[fID].GetSize())) { return false; }
      if (prog[fID].GetTags().size() != prog_func_num_tags) { return false; }
      for (size_t iID = 0; iID < prog[fID].GetSize(); ++iID) {
        if (prog[fID][iID].GetArgs().size() != prog_inst_num_args) { return false; }
        if (prog[fID][iID].GetTags().size() != prog_inst_num_tags) { return false; }
        for (size_t k = 0; k < prog[fID][iID].GetArgs().size(); ++k) {
          if (!prog_inst_arg_val_range.Valid(prog[fID][iID].GetArg(k))) { return false; }
        }
      }
    }
    return true;
  }
};

}