#pragma once

#include <utility>
#include <unordered_map>

#include "emp/base/Ptr.hpp"
#include "emp/datastructs/map_utils.hpp"
#include "emp/tools/string_utils.hpp"

#include "psb/readers/FizzBuzz.hpp"

#include "../BaseProblemHardware.hpp"
#include "../Event.hpp"
#include "../TestResult.hpp"
#include "BaseProblem.hpp"

namespace psynth::problems {

struct FizzBuzzHardware : public BaseProblemHardware {

  enum class CATEGORY { FIZZ, BUZZ, FIZZBUZZ, ECHO, NONE };

  CATEGORY out_category = CATEGORY::NONE;
  int echo_num = 0;

  void Reset() override {
    out_category = CATEGORY::NONE;
    echo_num = 0;
  }

  void SubmitFizz() {
    out_category = CATEGORY::FIZZ;
  }

  void SubmitBuzz() {
    out_category = CATEGORY::BUZZ;
  }

  void SubmitFizzBuzz() {
    out_category = CATEGORY::FIZZBUZZ;
  }

  void SubmitEcho(int value) {
    out_category = CATEGORY::ECHO;
    echo_num = value;
  }

  CATEGORY GetOutputCategory() {
    return out_category;
  }

  int GetOutputEchoNumber() {
    return echo_num;
  }

};

struct FizzBuzz : public BaseProblem {
  using reader_t = psb::readers::FizzBuzz;
  using input_t = typename reader_t::input_t;
  using output_t = typename reader_t::output_t;
  using test_case_t = std::pair<input_t, output_t>;
  using prob_hw_t = FizzBuzzHardware;

  size_t input_sig_event_id = 0;
  double max_test_score = 1.0;
  std::unordered_map<
    std::string,
    FizzBuzzHardware::CATEGORY
  > output_str_to_category = {
    {"Fizz", FizzBuzzHardware::CATEGORY::FIZZ},
    {"Buzz", FizzBuzzHardware::CATEGORY::BUZZ},
    {"FizzBuzz", FizzBuzzHardware::CATEGORY::FIZZBUZZ}
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
      "SubmitFizz",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitFizz();
      },
      "Categorize test input as Fizz"
    );

    inst_lib.AddInst(
      "SubmitBuzz",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitBuzz();
      },
      "Categorize test input as Buzz"
    );

    inst_lib.AddInst(
      "SubmitFizzBuzz",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        component.SubmitFizzBuzz();
      },
      "Categorize test input as FizzBuzz"
    );

    inst_lib.AddInst(
      "SubmitEcho",
      [](hardware_t& hw, const inst_t& inst) {
        auto& component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
        auto& call_state = hw.GetCurThread().GetExecState().GetTopCallState();
        auto& mem_state = call_state.GetMemory();
        int output = (int)mem_state.AccessWorking(inst.GetArg(0));
        component.SubmitEcho(output);
      },
      "Echo input as output"
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

    // Get correct output
    // auto correct_output = output_str_to_category[correct_output_str];
    FizzBuzzHardware::CATEGORY correct_output = (emp::Has(output_str_to_category, correct_output_str)) ?
      output_str_to_category[correct_output_str] :
      FizzBuzzHardware::CATEGORY::ECHO;

    // Get reference to the problem component
    auto& prob_component = hw.GetCustomComponent().template GetProbHW<prob_hw_t>();
    // Did the program record any output?
    const bool has_output = prob_component.GetOutputCategory() != FizzBuzzHardware::CATEGORY::NONE;
    // Is the output correct category?
    bool correct = prob_component.GetOutputCategory() == correct_output;
    // If output is correct and correct output is echo, check echo number.
    if (correct && (correct_output == FizzBuzzHardware::CATEGORY::ECHO)) {
      const int correct_echo = emp::from_string<int>(correct_output_str);
      correct = prob_component.GetOutputEchoNumber() == correct_echo;
    }
    emp_assert( !correct || (has_output && correct) ); // If it's correct, it must have output.
    // No partial credit on this problem.
    const double partial_credit = (correct) ? max_test_score : 0.0;
    return {has_output, correct, partial_credit};
  }

};

}