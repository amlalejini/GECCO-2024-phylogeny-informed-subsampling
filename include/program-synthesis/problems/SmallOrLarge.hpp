#pragma once

#include <utility>

#include "emp/base/Ptr.hpp"

#include "psb/readers/SmallOrLarge.hpp"

#include "../BaseProblemHardware.hpp"
#include "BaseProblem.hpp"
#include "../Event.hpp"

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

};

struct SmallOrLarge : public BaseProblem {
  using reader_t = psb::readers::SmallOrLarge;
  using input_t = typename reader_t::input_t;
  using output_t = typename reader_t::output_t;
  using test_case_t = std::pair<input_t, output_t>;
  using prob_hw_t = SmallOrLargeHardware;

  size_t input_sig_event_id = 0;

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
};

}