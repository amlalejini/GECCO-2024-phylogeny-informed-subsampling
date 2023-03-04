#pragma once

#include "BaseProblem.hpp"
#include "psb/readers/Median.hpp"

namespace psynth::problems {

struct Median : public BaseProblem {
  using reader_t = psb::readers::Median;

};

}