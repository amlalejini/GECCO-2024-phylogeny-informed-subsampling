#pragma once

namespace psynth {

struct BaseProblemHardware {

  virtual ~BaseProblemHardware() = default;
  virtual void Reset() = 0;

};

struct NumericOutputHardware : public BaseProblemHardware {
  double output_value = 0.0;
  bool output_set = false;

  void Reset() override {
    output_value = 0.0;
    output_set = false;
  }

  void SubmitOutput(double val) {
    output_set = true;
    output_value = val;
  }

  double GetOutput() const {
    return output_value;
  }

  bool HasOutput() const {
    return output_set;
  }

};

}