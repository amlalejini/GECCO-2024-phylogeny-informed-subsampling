#pragma once

#include "../BaseProblemHardware.hpp"
#include "BaseProblem.hpp"
#include "psb/readers/SmallOrLarge.hpp"

namespace psynth::problems {

struct SmallOrLargeHardware : public BaseProblemHardware {

  enum class CATEGORY { SMALL, LARGE, NEITHER, NONE };

  CATEGORY out_category = CATEGORY::NONE;

  void Reset() override {
    out_category = CATEGORY::NONE;
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
    // TODO
  }

  // Configure evaluation initialization

  // Configure

};

}