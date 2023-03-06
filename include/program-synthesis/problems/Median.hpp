#pragma once

#include "BaseProblem.hpp"
#include "psb/readers/Median.hpp"

namespace psynth::problems {

struct Median : public BaseProblem {
  using reader_t = psb::readers::Median;

  template<typename HARDWARE_T>
  void ConfigureHardware(HARDWARE_T& hw) { emp_assert(false); }

  template<typename INST_LIB_T>
  void AddInstructions(INST_LIB_T& inst_lib) { emp_assert(false); }

};

}