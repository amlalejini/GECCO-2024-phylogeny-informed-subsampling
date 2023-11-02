#pragma once

#include <utility>
#include <unordered_map>

#include "emp/base/Ptr.hpp"
#include "emp/datastructs/map_utils.hpp"

#include "psb/readers/Grade.hpp"

#include "../BaseProblemHardware.hpp"
#include "../Event.hpp"
#include "../TestResult.hpp"
#include "BaseProblem.hpp"

namespace psynth::problems {

struct GradeHardware : public BaseProblemHardware {

  enum class CATEGORY { A, B, C, D, F, NONE };

  CATEGORY out_category = CATEGORY::NONE;

  void Reset() override {
    out_category = CATEGORY::NONE;
  }

  void SubmitA() {
    out_category = CATEGORY::A;
  }

  void SubmitB() {
    out_category = CATEGORY::B;
  }

  void SubmitC() {
    out_category = CATEGORY::C;
  }

  void SubmitD() {
    out_category = CATEGORY::D;
  }

  void SubmitF() {
    out_category = CATEGORY::F;
  }

  CATEGORY GetOutput() {
    return out_category;
  }

};

struct Grade : public BaseProblem {
  using reader_t = psb::readers::Grade;
  using input_t = typename reader_t::input_t;
  using output_t = typename reader_t::output_t;
  using test_case_t = std::pair<input_t, output_t>;
  using prob_hw_t = GradeHardware;

  size_t input_sig_event_id = 0;
  double max_test_score = 1.0;
  std::unordered_map<
    std::string,
    GradeHardware::CATEGORY
  > output_str_to_category = {
    {"A", GradeHardware::CATEGORY::A},
    {"B", GradeHardware::CATEGORY::B},
    {"C", GradeHardware::CATEGORY::C},
    {"D", GradeHardware::CATEGORY::D},
    {"F", GradeHardware::CATEGORY::F}
  };

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
      "SubmitA",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitA();
      },
      "Categorize grade"
    );

    inst_lib.AddInst(
      "SubmitB",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitB();
      },
      "Categorize grade"
    );

    inst_lib.AddInst(
      "SubmitC",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitC();
      },
      "Categorize grade"
    );

    inst_lib.AddInst(
      "SubmitD",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitD();
      },
      "Categorize grade"
    );

    inst_lib.AddInst(
      "SubmitF",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitF();
      },
      "Categorize grade"
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
    const input_t& input = test_io.first;
    const size_t tag_w = HARDWARE_T::tag_t::GetSize();
    hw.QueueEvent(
      NumericMessageEvent<tag_w>{
        input_sig_event_id,
        hw.GetCustomComponent().GetInputTag(),
        {
          {0, (double)input[0]},
          {1, (double)input[1]},
          {2, (double)input[2]},
          {3, (double)input[3]},
          {4, (double)input[4]},
        }
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
    const bool has_output = prob_component.GetOutput() != GradeHardware::CATEGORY::NONE;
    // Is the output correct?
    const bool correct = prob_component.GetOutput() == correct_output;
    emp_assert( !correct || (has_output && correct) ); // If it's correct, it must have output.
    // No partial credit on this problem.
    const double partial_credit = (correct) ? max_test_score : 0.0;
    return {has_output, correct, partial_credit};
  }

};

}