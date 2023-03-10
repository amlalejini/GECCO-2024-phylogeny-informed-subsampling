#pragma once

#include <utility>
#include <unordered_map>

#include "emp/base/Ptr.hpp"
#include "emp/datastructs/map_utils.hpp"

#include "psb/readers/SmallOrLarge.hpp"

#include "../BaseProblemHardware.hpp"
#include "../Event.hpp"
#include "../TestResult.hpp"
#include "BaseProblem.hpp"

namespace psynth::problems {

struct SmallOrLargeHardware : public BaseProblemHardware {

  enum class CATEGORY { SMALL, LARGE, NEITHER, NONE };

  CATEGORY out_category = CATEGORY::NONE;

  void Reset() override {
    out_category = CATEGORY::NONE;
  }

  void SubmitSmall() {
    out_category = CATEGORY::SMALL;
  }

  void SubmitLarge() {
    out_category = CATEGORY::LARGE;
  }

  void SubmitNeither() {
    out_category = CATEGORY::NEITHER;
  }

  CATEGORY GetOutput() {
    return out_category;
  }

};

struct SmallOrLarge : public BaseProblem {
  using reader_t = psb::readers::SmallOrLarge;
  using input_t = typename reader_t::input_t;
  using output_t = typename reader_t::output_t;
  using test_case_t = std::pair<input_t, output_t>;
  using prob_hw_t = SmallOrLargeHardware;

  size_t input_sig_event_id = 0;
  std::unordered_map<
    std::string,
    SmallOrLargeHardware::CATEGORY
  > output_str_to_category = {
    {"small", SmallOrLargeHardware::CATEGORY::SMALL},
    {"", SmallOrLargeHardware::CATEGORY::NEITHER},
    {"large", SmallOrLargeHardware::CATEGORY::LARGE}
  };


  // Configure hardware
  template<typename HARDWARE_T>
  void ConfigureHardware(HARDWARE_T& hw) {
    auto& hw_component = hw.GetCustomComponent();
    hw_component.template CreateProblemHardware<prob_hw_t>();
    // TODO get event library, get name of requisite events
  }

  // Configure instructions
  template<typename INST_LIB_T>
  void AddInstructions(INST_LIB_T& inst_lib) {

    using hardware_t = typename INST_LIB_T::hardware_t;
    using inst_t = typename hardware_t::inst_t;

    inst_lib.AddInst(
      "SubmitSmall",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitSmall();
      },
      "Categorize test input as small"
    );

    inst_lib.AddInst(
      "SubmitLarge",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitLarge();
      },
      "Categorize test input as large"
    );

    inst_lib.AddInst(
      "SubmitNeither",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitNeither();
      },
      "Categorize test input as neither"
    );

  }

  template<typename EVENT_LIB_T>
  void AddEvents(EVENT_LIB_T& event_lib) {
    // Don't need to add events here, just record appropriate input signal event id
    input_sig_event_id = event_lib.GetID("NumericInputSignal");
  }

  // Configure evaluation initialization
  template<typename HARDWARE_T, typename ORG_T>
  void InitTest(HARDWARE_T& hw, ORG_T& org, const test_case_t& test_io) {
    // Load test inputs into hardware
    const input_t input = test_io.first;
    const size_t tag_w = HARDWARE_T::tag_t::GetSize();
    hw.QueueEvent(
      NumericMessageEvent<tag_w>{
        input_sig_event_id,
        hw.GetCustomComponent().GetInputTag(),
        {{0, (double)input}}
      }
    );
  }

  template<typename HARDWARE_T, typename ORG_T>
  psynth::TestResult EvaluateOutput(HARDWARE_T& hw, ORG_T& org, const test_case_t& test_io) {
    const output_t& correct_output_str = test_io.second;
    emp_assert(emp::Has(output_str_to_category, correct_output_str));
    // Get correct output
    auto correct_output = output_str_to_category[correct_output_str];
    // Get reference to the problem component
    auto& prob_component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
    // Did the program record any output?
    const bool has_output = prob_component.GetOutput() != SmallOrLargeHardware::CATEGORY::NONE;
    // Is the output correct?
    const bool correct = prob_component.GetOutput() == correct_output;
    emp_assert( !correct || (has_output && correct) ); // If it's correct, it must have output.
    // No partial credit on this problem.
    const double partial_credit = (correct) ? 1.0 : 0.0;
    return {has_output, correct, partial_credit};
  }

};

}