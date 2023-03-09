#pragma once

#include <utility>

#include "BaseProblem.hpp"
#include "psb/readers/Median.hpp"

namespace psynth::problems {

struct Median : public BaseProblem {
  using reader_t = psb::readers::Median;
  using input_t = typename reader_t::input_t;
  using output_t = typename reader_t::output_t;
  using test_case_t = std::pair<input_t, output_t>;

  template<typename HARDWARE_T>
  void ConfigureHardware(HARDWARE_T& hw) { emp_assert(false); }

  template<typename INST_LIB_T>
  void AddInstructions(INST_LIB_T& inst_lib) { emp_assert(false); }

  template<typename EVENT_LIB_T>
  void AddEvents(EVENT_LIB_T& event_lib) { emp_assert(false); }

  template<typename HARDWARE_T, typename ORG_T>
  void InitTest(HARDWARE_T& hw, ORG_T& org, const test_case_t& test_io) { emp_assert(false); }

};

}