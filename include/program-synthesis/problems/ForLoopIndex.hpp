#pragma once

#include <utility>
#include <unordered_map>

#include "emp/base/Ptr.hpp"
#include "emp/math/sequence_utils.hpp"

#include "psb/readers/ForLoopIndex.hpp"

#include "../BaseProblemHardware.hpp"
#include "../Event.hpp"
#include "../TestResult.hpp"
#include "BaseProblem.hpp"

namespace psynth::problems {

struct ForLoopIndexHardware : public BaseProblemHardware {
  emp::vector<int> output;

  void Reset() override {
    ClearOutput();
  }

  void SubmitOutput(int value) {
    output.emplace_back(value);
  }

  void ClearOutput() {
    output.clear();
  }

  bool HasOutput() const {
    return output.size() > 0;
  }

  const emp::vector<int>& GetOutput() const {
    return output;
  }
};

struct ForLoopIndex : public BaseProblem {
  using reader_t = psb::readers::ForLoopIndex;
  using input_t = typename reader_t::input_t;
  using output_t = typename reader_t::output_t;
  using test_case_t = std::pair<input_t, output_t>;
  using prob_hw_t = ForLoopIndexHardware;

  size_t input_sig_event_id = 0;
  double max_test_score = 1.0;

  // Configure hardware
  template<typename HARDWARE_T>
  void ConfigureHardware(HARDWARE_T& hw) {
    auto& hw_component = hw.GetCustomComponent();
    hw_component.template CreateProblemHardware<prob_hw_t>();
  }

  // Configure instructions
  template<typename INST_LIB_T>
  void AddInstructions(INST_LIB_T& inst_lib) {

    using hardware_t = typename INST_LIB_T::hardware_t;
    using inst_t = typename hardware_t::inst_t;

    inst_lib.AddInst(
      "SubmitOutput",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        auto& call_state = hw.GetCurThread().GetExecState().GetTopCallState();
        auto& mem_state = call_state.GetMemory();
        const double output = mem_state.AccessWorking(inst.GetArg(0));
        component.SubmitOutput(output);
      },
      "Print output"
    );

    inst_lib.AddInst(
      "ClearOutput",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.ClearOutput();
      },
      "Reset output buffer"
    );

  }

  template<typename EVENT_LIB_T>
  void AddEvents(EVENT_LIB_T& event_lib) {
    input_sig_event_id = event_lib.GetID("NumericInputSignal");
  }

  template<typename HARDWARE_T, typename ORG_T>
  void InitTest(HARDWARE_T& hw, ORG_T& org, const test_case_t& test_io) {
    // Load test inputs into hardware
    const input_t& input = test_io.first;
    const size_t tag_w = HARDWARE_T::tag_t::GetSize();
    hw.QueueEvent(
      NumericMessageEvent<tag_w>{
        input_sig_event_id,
        hw.GetCustomComponent().GetInputTag(),
        {
          {0, (double)input[0]},
          {1, (double)input[1]},
          {2, (double)input[2]}
        }
      }
    );
  }

  template<typename HARDWARE_T, typename ORG_T>
  TestResult EvaluateOutput(HARDWARE_T& hw, ORG_T& org, const test_case_t& test_io) {
    // Get correct output
    const auto& correct_output = test_io.second;
    // Get reference to problem component
    auto& prob_component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
    const auto& output = prob_component.GetOutput();
    // Did the program record any output?
    const bool has_output = prob_component.HasOutput();
    if (!has_output) {
      return {has_output, false, 0};
    }

    // Is the output correct?
    const bool correct = (output == correct_output);
    if (correct) {
      return {has_output, correct, max_test_score};
    }

    // Have output, but not correct.
    const double max_dist = emp::Max(correct_output.size(), output.size());
    const double dist = emp::calc_edit_distance(correct_output, output);
    emp_assert(max_dist > 0);
    const double modifier = (max_dist - dist) / max_dist;
    return {has_output, correct, max_test_score * modifier};
  }


};

}