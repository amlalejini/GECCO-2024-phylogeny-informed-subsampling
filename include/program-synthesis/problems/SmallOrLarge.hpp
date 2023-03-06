#pragma once

#include <utility>

#include "emp/base/Ptr.hpp"

#include "psb/readers/SmallOrLarge.hpp"

#include "../BaseProblemHardware.hpp"
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

};

struct SmallOrLarge : public BaseProblem {
  using reader_t = psb::readers::SmallOrLarge;
  using input_t = typename reader_t::input_t;
  using output_t = typename reader_t::output_t;
  using prob_hw_t = SmallOrLargeHardware;

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

  // Configure evaluation initialization

  // Configure

};

}