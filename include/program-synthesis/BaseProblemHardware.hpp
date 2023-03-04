#pragma once

namespace psynth {

struct BaseProblemHardware {

  virtual ~BaseProblemHardware() = default;
  virtual void Reset() = 0;

};

}